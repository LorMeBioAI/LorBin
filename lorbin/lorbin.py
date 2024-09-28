import argparse
import logging
import os
import sys
from Bio import SeqIO
from .generate_coverage import generate_cov,combine_cov
from .generate_kmer import generate_kmer_features_from_fasta
from .fasta import fasta_iter
from multiprocessing.pool import Pool
from .atomicwrite import atomic_write
import logging
import pandas as pd
from . import utils as utils
from .model import vae as vae
import torch
import numpy as np
from torch.utils.data.dataset import TensorDataset
from torch.utils.data import DataLoader
from .cluster import bin_cluster
from .check_arg import check_generate_data, check_cluster, check_train


def concat(outdir, fastas):
    fw = open(f'{outdir}','w')
    for e,fasta in enumerate(fastas):
        for record in SeqIO.parse(fasta, 'fasta'):
            fw.write(f'>{e}-{record.id}\n')
            fw.write(f'{record.seq}\n')
    fw.close()

def generate_data(logger,fastadir, bams, outdir, num_process):
    # generate kmer and abundance data
    if not check_generate_data(logger, fastadir, bams, outdir):
        sys.exit()
    logger.info('process kmer')
    kmer, length = generate_kmer_features_from_fasta(fastadir)
    length.to_csv(f'{outdir}/length.csv')
    logger.info('process bams')
    with Pool(num_process if num_process != 0 else None) as pool:
        results = [
             pool.apply_async(
                generate_cov,
                args=(
                    bam_file,
                    bam_index,
                    outdir,
                    logger
                ))
            for bam_index, bam_file in enumerate(bams)]
        for r in results:
            s = r.get()
    data_cov = combine_cov(outdir, bams)
    data = pd.merge(kmer, data_cov, how='inner', on=None,left_index=True, right_index=True, sort=False, copy=True)
    data = pd.merge(data, length, how='inner', on=None,left_index=True, right_index=True, sort=False, copy=True)
    with atomic_write(os.path.join(outdir, 'data.csv'), overwrite=True) as ofile:
        data.to_csv(ofile)

def generate_markers(logger,fastadir, bin_length, num_process, outdir):
    logger.info('to get single-copy markers')
    utils.generate_markers(fastadir, bin_length, num_process, output = outdir)


def train_vae(logger, outdir, batch_size=64, epoch=300, batchsteps=[25],  lrate=1e-3, datapath=None):
    if not datapath:
        datapath = f"{outdir}/data.csv"
    df = pd.read_csv(datapath, index_col=0)
    data = df.values
    tnf = np.array(data[:,:136],dtype=np.float32)
    rpkm = np.array(data[:,136:-1],dtype=np.float32)
    length = data[:,-1]
    rpkm, tnf, length, weights = vae.normalize(rpkm, tnf,length)
    dataloader = vae.make_dataloader(rpkm, tnf, length, weights, batch_size)
    is_cuda = torch.cuda.is_available()
    logger.info(f"VAE build cuda set is {is_cuda}")
    model = vae.VAE(rpkm.shape[1], None, 32, None, 200.0, 0.2, is_cuda,0)
    model.trainmodel(logger, outdir,dataloader, nepochs=epoch, batchsteps = batchsteps)
    testdataloader = DataLoader(dataset=dataloader.dataset,batch_size=batch_size,drop_last=False,shuffle=False)
    latent = model.get_latent(testdataloader)
    pd.DataFrame(latent,index=df.index).to_csv(f'{outdir}/embedding.csv')

def cluster(logger, outdir,fastadir,embeddingdir, bin_length, feature,a):
    df = pd.read_csv(embeddingdir,index_col=0)
    names = df.index
    contig_all = names
    contig_list = list(names)
    contig2marker = utils.get_marker(f"{outdir}/markers.hmmout",orf_finder="prodigal",min_contig_len=1500, contig_names= contig_all, fasta_path=fastadir)
    if not os.path.exists(embeddingdir):
        logger.info("No embedding!")
        sys.exit()
    latent = df.values
    if not contig2marker:
        logger.info("no markers")
    contig_dict = utils.process_fasta(fastadir)
    labels, keep = bin_cluster(logger, latent, contig2marker, contig_dict, contig_list, contig_all, bin_length, feature, a)
    pd.DataFrame({'label':labels},index=contig_all).to_csv(f'{outdir}/label.csv')
    write_bin(contig_all, labels,contig_dict,f"{outdir}/output_bins",bin_length)

def mcluster(logger, outdir, fastadir, embeddingdir, bin_length, feature,a):
    if generate_cluster(logger, outdir, fastadir, embeddingdir):
        sys.exit()
    df = pd.read_csv(embeddingdir,index_col=0)
    names = df.index
    contig_all = names
    contig_list = list(names)
    nsample, samplenames = msample(names)
    contig2marker = utils.get_marker(f"{outdir}/markers.hmmout",orf_finder="prodigal",min_contig_len=1500, contig_names= contig_all, fasta_path=fastadir)
    if not os.path.exists(embeddingdir):
        logger.info("No embedding!")
        sys.exit()
    if not contig2marker:
        logger.info("no markers")
    contig_dict = utils.process_fasta(fastadir)
    embedding = df.values
    for i in range(len(nsample)-1):
        latent = embedding[nsample[i]:nsample[i+1]]
        labels, keep = bin_cluster(logger, latent, contig2marker, contig_dict, contig_list, contig_all, bin_length, feature, a)
        pd.DataFrame({'label':labels},index=contig_all).to_csv(f'{outdir}/label_{samplename[i]}.csv')
        write_bin(contig_all, labels,contig_dict,f"{output}/output_bins_{samplename[i]}",bin_length)


def msample(names):
    samplestart=[]
    samplename=[]
    prefix=" "
    for e,name in enumerate(names):
        sample = name.split('-')[0]
        if sample!=prefix:
            samplestart.append(e)
            prefix=sample
            samplename.sppend(prefix)
    samplestart.append(e)
    return samplestart, samplename

def write_bin(namelist,contig_labels,contig_seqs,output, minfasta):
    from collections import defaultdict

    res = defaultdict(list)
    for label, name in zip(contig_labels, namelist):
        if label != -1:
            res[label].append(name)

    os.makedirs(output, exist_ok=True)

    for label, contigs in res.items():
        sizes = [len(contig_seqs[contig]) for contig in contigs]
        whole_bin_bp = sum(sizes)
        if whole_bin_bp >= minfasta:
            ofname = f'bin.{label}.fa'
            ofname = os.path.join(output, ofname)
            with open(ofname,'w') as ofile:
                for contig in contigs:
                    ofile.write(f'>{contig}\n{contig_seqs[contig]}\n')

def parser_args():
    parser = argparse.ArgumentParser(
        #description=doc,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=True,)
    subparsers = parser.add_subparsers(title="LorBin subcommands", dest='cmd', metavar='')
    bin_mode = subparsers.add_parser('bin',help="Bin contigs using one command.")
    generate_data = subparsers.add_parser('generate_data', help="Generate sequence features(kmer and abundance) as trainning data.")
    concat_fasta = subparsers.add_parser('concat', help="""Create the input FASTA file for LorBin, input should be at least two FASTA files, each from a sample-specific assembly, resulting FASTA can be binsplit with separator '-'""")
    train = subparsers.add_parser('train', help="Train model and get embedding.")
    cluster = subparsers.add_parser('cluster', help='Bin contigs using two-stage clustering algorithm based on embedded features file.')
    for p in [ bin_mode,generate_data,cluster]:
        p.add_argument('-o','--output',type=str,help='Output directory (will be created if non-existent)',required=True, default=None)
        p.add_argument('-fa','--fasta',type=str, help='Path to the input fasta file.',required=True, default=None)
        p.add_argument('--bin_length', default=80000,help='Minimum bin size in bps (Default: 80000)')
    for p in [ bin_mode, generate_data]:
        p.add_argument('-b','--bam',type=str, nargs='+',help='Path to the input BAM(.bam) file. ',required=True,default=None)
        p.add_argument('--num_process', default=10, help=' Number of threads used (default: 10)')
    for p in [bin_mode, cluster]:
        p.add_argument('--evaluation',type=str, default="no_markers", help='Evaluation model used(no_markers, markers110, markers35, default: nomarkers')
        p.add_argument('-a','--akeep',default=0.6, help='The cut-off parameters of re-clustering decision model(0~1, default:0.6)')
        p.add_argument('--multi',action='store_true', default=False, help='Cluster uses more samples')
    cluster.add_argument('--embeddingdir','-e',required=True,default=None, help='The path of embedding csv file used in clustering')
    cluster.add_argument('--num_process', default=10, help=' Number of threads used (default: 10)')
    train.add_argument('--data',type=str, help='The path of training data', required=True)
    train.add_argument('-o','--output',type=str,help='Output directory (will be created if non-existent)',required=True, default=None)
    train.add_argument('--epoch','-n', help='training epoch (default: 300)', default=50)
    train.add_argument('--lrate','-l',help='learning rate (default: 0.001)', default=0.001)
    train.add_argument('--batch_size', help = 'batch size (default: 64)', default=64)
    train.add_argument('--batchsteps', help = 'batchseteps (default: 30, 60, 120)', default=[10, 20], nargs='+')
    concat_fasta.add_argument('-fa','--fasta',type=str, nargs='+',help='The path to input FASTA files',required=True)
    concat_fasta.add_argument('-o','--output',help="The path to output FASTA file", required=True)
    args = parser.parse_args()
    return args

def main():
    args=parser_args()
    logger = logging.getLogger('LorBin')
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")# Set the log output format
    if args.cmd=="bin":
        if not check_generate_data(logger, args.fasta, args.bam, args.output):sys.exit()
        file_handler = logging.FileHandler(f"{args.output}/LorBin.log") # Create a Handler object to control where the log is output
        file_handler.setFormatter(formatter) 
        logger.addHandler(file_handler)
        generate_data(logger, args.fasta, args.bam, args.output, args.num_process)
        generate_markers(logger, args.fasta, args.bin_length, args.num_process, args.output)
        train_vae(logger,args.output, epoch=30)
        embeddingdir = f"{args.output}/embedding.csv"
        if args.multi:
            mcluster(logger, args.output, args.fasta, embeddingdir, args.bin_length, args.evaluation,args.akeep)
        else:
            cluster(logger, args.output, args.fasta, embeddingdir, args.bin_length, args.evaluation, args.akeep)
    elif args.cmd=='generate_data':
        if not check_generate_data(logger, args.fasta, args.bam, args.output):sys.exit()
        file_handler = logging.FileHandler(f"{args.output}/LorBin.log") 
        logger.addHandler(file_handler)
        file_handler.setFormatter(formatter)
        generate_data(logger, args.fasta, args.bam, args.output, args.num_process)
        #generate_markers(logger, args.fasta, args.bin_length, args.num_process, args.output)
    elif args.cmd == 'train':
        check_train(logger, args.output)
        file_handler = logging.FileHandler(f"{args.output}/LorBin.log") 
        logger.addHandler(file_handler)
        file_handler.setFormatter(formatter) 
        train_vae(logger, args.output,  args.batch_size, args.epoch, args.batchsteps, args.lrate, args.data)
    elif args.cmd=='cluster':
        if not check_cluster(logger, args.output, args.fasta, args.embeddingdir, args.evaluation, args.akeep):sys.exit()
        file_handler = logging.FileHandler(f"{args.output}/LorBin.log") 
        logger.addHandler(file_handler)
        file_handler.setFormatter(formatter)
        generate_markers(logger, args.fasta, args.bin_length, args.num_process, args.output)
        if args.multi:
            mcluster(logger, args.output, args.fasta, args.embeddingdir, args.bin_length, args.evaluation,args.akeep)
        else:
            cluster(logger, args.output, args.fasta, args.embeddingdir, args.bin_length, args.evaluation,args.akeep)
    elif args.cmd=='concat':
        concat(args.output,args.fasta)
    else:
        logger.info('cmd error!')
        sys.exit()
