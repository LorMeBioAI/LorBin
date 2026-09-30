from collections import defaultdict
import pandas as pd
import torch
from Bio import SeqIO
from sklearn.neighbors import kneighbors_graph
import math
from sklearn.feature_extraction.text import CountVectorizer, TfidfTransformer
from .utils import get_marker, process_fasta
from .model.KeepModel import KeepModel
from .model.EvaluationModel import EvaluationModel
import numpy as np
from sklearn.cluster import Birch, DBSCAN
import logging
import os
import time



class CandidatePool:
    """The candidate bins that get_bin_best scans, with each candidate's length,
    marker count and distinct-marker count kept up to date as contigs are
    extracted, so a pass does not re-walk every candidate's contigs.

    best() and extract() pick and remove exactly what get_bin_best and the
    removal loops in bin_cluster do, in the same pool order: F1 depends only on
    (distinct markers, total markers), so each distinct pair is scored once
    with the same single-row model call get_bin_best makes.
    """

    def __init__(self, resultpool, contig_len, contig_markers):
        n_cand = len(resultpool)
        sizes = np.fromiter((len(r) for r in resultpool), dtype=np.int64, count=n_cand)
        self.ptr = np.zeros(n_cand + 1, dtype=np.int64)
        np.cumsum(sizes, out=self.ptr[1:])
        self.members = np.fromiter((c for r in resultpool for c in r), dtype=np.int64,
                                   count=int(self.ptr[-1]))
        owner = np.repeat(np.arange(n_cand, dtype=np.int64), sizes)
        # contig -> candidates containing it
        order = np.argsort(self.members, kind='stable')
        self.inv_cand = owner[order]
        self.inv_ptr = np.searchsorted(self.members[order], np.arange(len(contig_len) + 1))
        self.alive = np.ones(len(contig_len), dtype=bool)
        self.contig_len = contig_len

        cum = np.zeros(len(self.members) + 1, dtype=np.int64)
        np.cumsum(contig_len[self.members], out=cum[1:])
        self.weight = cum[self.ptr[1:]] - cum[self.ptr[:-1]]

        # per-candidate marker counts, only for candidates holding a marker contig
        self.contig_markers = contig_markers
        marker_contigs = [c for c in contig_markers if self.inv_ptr[c + 1] > self.inv_ptr[c]]
        self.row = np.full(n_cand, -1, dtype=np.int64)
        if marker_contigs:
            has = np.unique(np.concatenate([self._cands(c) for c in marker_contigs]))
            self.row[has] = np.arange(len(has))
        n_markers = 1 + max((m for ms in contig_markers.values() for m in ms), default=0)
        self.counts = np.zeros((int((self.row >= 0).sum()), n_markers), dtype=np.int32)
        for c in marker_contigs:
            rows = self.row[self._cands(c)]
            for m in contig_markers[c]:
                self.counts[rows, m] += 1
        self.n = np.zeros(n_cand, dtype=np.int64)
        self.u = np.zeros(n_cand, dtype=np.int64)
        self.f1 = np.zeros(n_cand, dtype=np.float64)
        self.contamination = np.zeros(n_cand, dtype=np.float64)
        self._scores = {}

    def _cands(self, c):
        return self.inv_cand[self.inv_ptr[c]:self.inv_ptr[c + 1]]

    def score_all(self, evaluationmodel):
        self.evaluationmodel = evaluationmodel
        cands = np.flatnonzero(self.row >= 0)
        self._rescore(cands)

    def _rescore(self, cands):
        rows = self.counts[self.row[cands]]
        self.n[cands] = rows.sum(axis=1)
        self.u[cands] = (rows > 0).sum(axis=1)
        for k in cands:
            n, u = int(self.n[k]), int(self.u[k])
            if n == 0:
                self.f1[k] = 0.0
                self.contamination[k] = 0.0
                continue
            # the same float arithmetic as get_bin_best
            contamination = (n - u) / n
            self.contamination[k] = contamination
            f1 = self._scores.get((u, n))
            if f1 is None:
                recall = u / 107
                f1 = self.evaluationmodel(torch.Tensor([[recall, 1 - contamination, abs(recall - 1 + contamination)]])).item()
                self._scores[(u, n)] = f1
            self.f1[k] = f1

    def best(self, keepmodel, minfasta, a):
        """get_bin_best: (candidate index, keep) or None."""
        eligible = np.flatnonzero((self.weight >= minfasta) & (self.n > 0))
        for max_contamination in [0.1, 0.2, 0.3, 0.4, 0.5, 1]:
            sel = eligible[self.contamination[eligible] <= max_contamination]
            if sel.size == 0:
                continue
            f1 = self.f1[sel]
            max_f1 = f1.max()
            if max_f1 <= 0:
                continue
            # highest F1; among ties the lightest, and among those the last in pool order
            tied = sel[f1 == max_f1]
            w = self.weight[tied]
            k = int(tied[np.flatnonzero(w == w.min())[-1]])
            n, u = int(self.n[k]), int(self.u[k])
            recall = u / 107
            contamination = (n - u) / n
            keep = keepmodel(torch.Tensor([[recall, 1 - contamination, abs(recall - 1 + contamination)]]))
            return k, (True if keep > a else False)
        return None

    def contigs(self, k):
        m = self.members[self.ptr[k]:self.ptr[k + 1]]
        return m[self.alive[m]].tolist()

    def extract(self, contigs):
        """Remove these contigs from every candidate."""
        contigs = np.asarray(contigs, dtype=np.int64)
        self.alive[contigs] = False
        touched_markers = []
        for c in contigs.tolist():
            cands = self._cands(c)
            self.weight[cands] -= self.contig_len[c]
            ms = self.contig_markers.get(c)
            if ms:
                rows = self.row[cands]
                for m in ms:
                    self.counts[rows, m] -= 1
                touched_markers.append(cands)
        if touched_markers:
            self._rescore(np.unique(np.concatenate(touched_markers)))

    def eligible_count(self, minfasta):
        return int(((self.weight >= minfasta) & (self.n > 0)).sum())


def _pool_inputs(contig2marker, namelist, contig_dict):
    contig_len = np.array([len(contig_dict[name]) for name in namelist], dtype=np.int64)
    marker_ids = {}
    contig_markers = {}
    for i, name in enumerate(namelist):
        ms = contig2marker.get(name)
        if ms:
            contig_markers[i] = [marker_ids.setdefault(m, len(marker_ids)) for m in ms]
    return contig_len, contig_markers


def extract_bins(logger, pool, keepmodel, evaluationmodel, minfasta, a, max_rejects, stage):
    """The extraction loop of bin_cluster, on a CandidatePool.

    Returns (extracted, all_extracted, keep_label): kept bins, every chosen bin,
    and the keep decision for each chosen bin. max_rejects is the run of
    consecutive rejected bins that stops the loop, or None for no limit.
    """
    pool.score_all(evaluationmodel)
    extracted, all_extracted, keep_label = [], [], []
    rejects = 0
    last_log = time.time()
    while True:
        found = pool.best(keepmodel, minfasta, a)
        if found is None or (max_rejects is not None and rejects > max_rejects):
            break
        k, keep = found
        max_bin = pool.contigs(k)
        if keep or max_rejects is None:
            extracted.append(max_bin)
        if keep:
            rejects = 0
        else:
            rejects += 1
        all_extracted.append(max_bin)
        keep_label.append(keep)
        pool.extract(max_bin)
        if time.time() - last_log > 600:
            logger.info(f"{stage}: {len(all_extracted)} bins extracted, "
                        f"{pool.eligible_count(minfasta)} candidates still eligible")
            last_log = time.time()
    logger.info(f"{stage}: {len(all_extracted)} bins extracted")
    return extracted, all_extracted, keep_label


def get_bin_best(keepmodel,evaluationmodel, resultpool, contig_to_marker, namelist, contig_dict, minfasta,a=0.6):
    for max_contamination in [0.1, 0.2, 0.3,0.4, 0.5,1]:
        max_F1 = 0
        weight_of_max = 1e9
        max_bin = None
        bin_recall = 0
        bin_contamination = 0

        for bin_contig_index in resultpool:
            cur_weight = sum(len(contig_dict[namelist[contig_index]]) for contig_index in bin_contig_index)
            if cur_weight < minfasta:
                continue
            marker_list = []
            for contig_index in bin_contig_index:
                marker_list.extend(contig_to_marker[namelist[contig_index]])
            if len(marker_list) == 0:
                continue
            recall = len(set(marker_list)) / 107
            contamination = (len(marker_list) - len(set(marker_list))) / len(marker_list)
            if contamination <= max_contamination:
                F1 = (evaluationmodel(torch.Tensor([[recall, 1 - contamination, abs(recall - 1 + contamination)]])))
                if F1 > max_F1:
                    max_F1 = F1
                    weight_of_max = cur_weight
                    max_bin = bin_contig_index
                    bin_recall = recall
                    bin_contamination = contamination
                elif F1 == max_F1 and cur_weight <= weight_of_max:
                    weight_of_max = cur_weight
                    max_bin = bin_contig_index
                    bin_recall = recall
                    bin_contamination = contamination
        if max_F1 > 0:  # if there is a bin with F1 > 0
            keep = keepmodel(torch.Tensor([[bin_recall, 1 - bin_contamination, abs(bin_recall - 1 + bin_contamination)]]))
            keep = True if keep > a else False
            return max_bin, keep

def get_bin_best_markers(keepmodel,evaluationmodel, vectorizer, tfidf_transformer, resultpool, contig_to_marker, namelist, contig_dict, minfasta,a=0.6):
    for max_contamination in [0.1, 0.2, 0.3,0.4, 0.5,1]:
        max_F1 = 0
        weight_of_max = 1e9
        max_bin = None
        bin_recall = 0
        bin_contamination = 0

        for bin_contig_index in resultpool:
            cur_weight = sum(len(contig_dict[namelist[contig_index]]) for contig_index in bin_contig_index)
            if cur_weight < minfasta:
                continue
            marker_list = []
            for contig_index in bin_contig_index:
                marker_list.extend(contig_to_marker[namelist[contig_index]])
            if len(marker_list) == 0:
                continue
            sentence = ' '.join(marker_list)
            X = vectorizer.transform([sentence])
            X = X.toarray()
            X_tfidf = tfidf_transformer.fit_transform(X)
            X = X_tfidf.toarray()

            recall = len(set(marker_list)) / 107
            contamination = (len(marker_list) - len(set(marker_list))) / len(marker_list)
            if contamination <= max_contamination:
                x_in1 = np.array([[recall, 1 - contamination, abs(recall - 1 + contamination)]])
                x_in = np.concatenate((x_in1, X), axis=1)
                x_in = torch.from_numpy(x_in).float()
                F1 = evaluationmodel(x_in)
                if F1 > max_F1:
                    max_F1 = F1
                    weight_of_max = cur_weight
                    max_bin = bin_contig_index
                    bin_recall = recall
                    bin_contamination = contamination
                    bin_x = X
                elif F1 == max_F1 and cur_weight <= weight_of_max:
                    weight_of_max = cur_weight
                    max_bin = bin_contig_index
                    bin_recall = recall
                    bin_contamination = contamination
                    bin_x = X
        if max_F1 > 0:
            x_in1 = np.array([[bin_recall, 1 - bin_contamination, abs(bin_recall - 1 + bin_contamination)]])
            x_in = np.concatenate((x_in1, bin_x), axis=1)
            x_in = torch.from_numpy(x_in).float()
            keep = keepmodel(x_in)
            keep = True if keep > a else False
            return max_bin, keep



def bin_cluster(logger, latent, contig2marker, contig_dict, contig_list, contig_all, minfasta, feature="no_markers", a=0.6):
    result_dict = {}
    # create BIRCH
    min_k_1 = min(200, latent.shape[0] - 1)
    dist_matrix = kneighbors_graph(
        latent,
        n_neighbors=min_k_1,
        mode='distance',
        p=2,
        n_jobs=10)
    dist_matrix_cos = kneighbors_graph(
        latent,
        n_neighbors=min_k_1,
        mode='distance',
        metric="cosine",
        p=2,
        n_jobs=10)


    cos_distance_data = dist_matrix_cos.data
    p2_distance_data = dist_matrix.data

    length_weight = np.array(
        [len(contig_dict[name]) for name in contig_list])
    resultpool = []
    contig_list_index = [i for i in range(len(contig_list))]

    eps_cos=[]
    eps_p2=[]
    for k in range(1,min(min_k_1,80)):
        eps_value=np.partition(cos_distance_data,latent.shape[0]*k-1)[latent.shape[0]*k-1]
        eps_cos.append(round(0.7*eps_value,2))
        eps_cos.append(round(1.2 * eps_value,2))
    eps_cos = list(set(eps_cos))
    eps_cos.sort()
    eps_cos = [ i for i in eps_cos if i<0.1] +[0.12,0.15,0.2]
    for k in range(1,min(80,min_k_1)):
        eps_value=np.partition(p2_distance_data,latent.shape[0]*k-1)[latent.shape[0]*k-1]
        eps_p2.append(math.floor(0.5*eps_value))
        eps_p2.append(math.floor(eps_value))
        eps_p2.append(math.floor(1.2 * eps_value))

    eps_p2 = list(set(eps_p2))
    eps_p2.sort()
    eps_p2 = eps_p2+[0.3*eps_p2[0],0.5*eps_p2[0],eps_p2[0]+0.5,eps_p2[1]+0.5]
    eps_p2.sort()


    for e,eps_value in enumerate(eps_p2):
        if eps_value<0.001:
            continue
        if e<5:
            minsample=5
        else:
            minsample=5
        dbscan = DBSCAN(eps=eps_value, min_samples=minsample, n_jobs=8, metric='precomputed')
        dbscan.fit(dist_matrix, sample_weight=length_weight)
        # predict the labal
        labels = dbscan.labels_
        res_temp = defaultdict(list)
        for label, name in zip(labels, contig_list_index):
            if label != -1:
                res_temp[label].append(name)
        resultpool.extend(res_temp.values())


    for e,eps_value in enumerate(eps_cos):
        if eps_value<0.001:
            continue
        if e<5:
            minsample=5
        else:
            minsample=5
        dbscan = DBSCAN(eps=eps_value, min_samples=minsample, n_jobs=8, metric='precomputed')
        dbscan.fit(dist_matrix_cos, sample_weight=length_weight)
        labels = dbscan.labels_
        result_dict[eps_value] = labels.tolist()

        res_temp = defaultdict(list)
        for label, name in zip(labels, contig_list_index):
            if label != -1:
                res_temp[f"cos_{label}"].append(name)
        resultpool.extend(res_temp.values())
    unique_resultpool = set(map(tuple, resultpool))
    resultpool = list(map(list, unique_resultpool))
    modeldir = os.path.split(__file__)[0]
    extracted = []
    all_extracted=[]
    if feature=="no_markers":
        model1 = EvaluationModel(3)
        model1.load_state_dict(torch.load(f'{modeldir}/model/huigui_weights.pt'))
        model2 = KeepModel(3)
        model2.load_state_dict(torch.load(f'{modeldir}/model/keepmodel_weights.pt'))
    elif feature=="markers110":
        model1 = EvaluationModel(110)
        model1.load_state_dict(torch.load(f'{modeldir}/model/huigui_weights_markers110.pt'))
        model2 = KeepModel(110)
        model2.load_state_dict(torch.load(f'{modeldir}/model/keepmodel_weights_markers110.pt'))
        vocabulary_df = pd.read_csv(f'{modeldir}/model/single-copy-genes.csv', index_col=0)
        vocabulary = list(vocabulary_df.iloc[:, 0])
        vectorizer = CountVectorizer(vocabulary=vocabulary, lowercase=False)
        tfidf_transformer = TfidfTransformer()
    else:
        model1 = EvaluationModel(35)
        model1.load_state_dict(torch.load(f'{modeldir}/model/huigui_weights_markers35.pt'))
        model2 = KeepModel(35)
        model2.load_state_dict(torch.load(f'{modeldir}/model/keepmodel_weights_markers35.pt'))
        vocabulary_df = pd.read_csv(f'{modeldir}/model/single-copy-genes35.csv', index_col=0)
        vocabulary = list(vocabulary_df.iloc[:, 0])
        vectorizer = CountVectorizer(vocabulary=vocabulary, lowercase=False)
        tfidf_transformer = TfidfTransformer()
    logger.info("load cluster quality asssement model and clustering decision model")
    keep_label=[]
    keep_count=0
    logger.info("cluster")
    total_length = sum(len(contig_dict[contig]) for contig in contig_list)
    if feature == "no_markers":
        contig_len, contig_markers = _pool_inputs(contig2marker, contig_all, contig_dict)
    fast = feature == "no_markers" and total_length >= minfasta and len(contig_list) != 1
    if fast:
        extracted, all_extracted, keep_label = extract_bins(
            logger, CandidatePool(resultpool, contig_len, contig_markers),
            model2, model1, minfasta, a, 200, "cluster")
    while not fast and total_length >= minfasta:
        if len(contig_list) == 1:
            extracted.append(contig_list_index)
            break
        if feature=="no_markers":
            max_bin = get_bin_best(model2, model1, resultpool, contig2marker, contig_all, contig_dict, minfasta,a)
        else:
            max_bin = get_bin_best_markers(model2, model1, vectorizer, tfidf_transformer,resultpool, contig2marker, contig_all, contig_dict, minfasta,a)
        if not max_bin or keep_count>200:
            break
        max_bin, keep = max_bin
        if keep:
            keep_count=0
            extracted.append(max_bin.copy())
        else:
            keep_count+=1
        all_extracted.append(max_bin.copy())
        keep_label.append(keep)
        for temp in max_bin.copy():
            for i, result in enumerate(resultpool):
                while temp in result:
                    result.remove(temp)

    contig2ix = {}
    label_index=0
    for i, cs in enumerate(extracted):
        for c in cs:
            contig2ix[contig_all[c]] = i
        label_index=i+1
    contig_labels = [contig2ix.get(c, -1) for c in contig_all]
    contig2ix_={}
    repeat=[]
    for i, cs in enumerate(all_extracted):
        for c in cs:
            if contig2ix_.get(contig_all[c], -1)==-1:
                contig2ix_[contig_all[c]] = i
            else:
                repeat.append(c)
        label_index=i+1
    contig_labels_ = [contig2ix_.get(c, -1) for c in contig_all]
    logger.info('start recluster')
    recluster_index = [index for index, value in enumerate(contig_labels) if value == -1]
    recluster_latent = latent[recluster_index]
    recluster_list = contig_all[recluster_index].tolist()
    min_k_2 = min(200, recluster_latent.shape[0]-1)

    dist_matrix = kneighbors_graph(
        recluster_latent,
        n_neighbors=min_k_2,
        mode='distance',
        p=2,
        n_jobs=10)



    p2_distance = dist_matrix.data
    eps_p2_2=[]


    resultpool = []
    for k in range(5,min(80,min_k_2),10):
        eps_value=np.partition(p2_distance,recluster_latent.shape[0]*k-1)[recluster_latent.shape[0]*k-1]
        eps_p2_2.append(math.floor(0.8*eps_value))
        eps_p2_2.append(math.floor(eps_value))
        eps_p2_2.append(math.floor(1.2*eps_value))
    eps_p2_2 = list(set(eps_p2_2))
    eps_p2_2.sort()

    threds =eps_p2_2 +[0.1*eps_p2_2[0],0.3*eps_p2_2[0],0.5*eps_p2_2[0],10,15,30]
    threds.sort()
    for thred in threds:
        if thred<0.00001:
            continue
        birch = Birch(threshold=thred, n_clusters=None)
        labels = birch.fit_predict(recluster_latent)

        res_temp = defaultdict(list)
        for label, name_index in zip(labels, recluster_index):
            if label != -1:
                res_temp[label].append(name_index)
        for res in res_temp.values():
            resultpool.append(res)
    unique_resultpool = set(map(tuple, resultpool))
    resultpool = list(map(list, unique_resultpool))

    minfasta = 1500
    keep_count=0
    recluster_length = sum(len(contig_dict[contig]) for contig in recluster_list)
    fast = feature == "no_markers" and recluster_length >= minfasta and len(recluster_list) != 1
    if fast:
        more, _, _ = extract_bins(logger, CandidatePool(resultpool, contig_len, contig_markers),
                                  model2, model1, minfasta, a, None, "recluster")
        extracted.extend(more)
    while not fast and recluster_length >= minfasta:
        if len(recluster_list) == 1:
            extracted.append(recluster_index)
            break

        if feature=="no_markers":
            max_bin = get_bin_best(model2, model1, resultpool, contig2marker, contig_all, contig_dict, minfasta,a)
        else:
            max_bin = get_bin_best_markers(model2, model1, vectorizer, tfidf_transformer,resultpool, contig2marker, contig_all, contig_dict, minfasta,a)
        if not max_bin or keep_count>100:
            break
        max_bin, keep = max_bin

        extracted.append(max_bin.copy())
        for temp in max_bin.copy():
            for i in range(len(resultpool)):
                if temp in set(resultpool[i]):
                    resultpool[i].remove(temp)
    contig2ix = {}
    for i, cs in enumerate(extracted):
        for c in cs:
            contig2ix[contig_all[c]] = i

    contig_labels = [contig2ix.get(c, -1) for c in contig_all]
    return contig_labels,keep_label

