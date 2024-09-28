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
    while sum(len(contig_dict[contig]) for contig in contig_list) >= minfasta:
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
    while sum(len(contig_dict[contig]) for contig in recluster_list) >= minfasta:
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

