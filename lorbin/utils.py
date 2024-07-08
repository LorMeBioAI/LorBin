from .fasta import fasta_iter
import subprocess
import multiprocessing
import sys
import contextlib
import tempfile
import os

normalize_marker_trans__dict = {
    'TIGR00388': 'TIGR00389',
    'TIGR00471': 'TIGR00472',
    'TIGR00408': 'TIGR00409',
    'TIGR02386': 'TIGR02387',
}

def get_marker(hmmout, contig_names, fasta_path=None, min_contig_len=None, orf_finder ="prodigal"):
    '''Parse HMM output file and return markers
    '''
    import pandas as pd
    data = pd.read_table(hmmout, sep=r'\s+',  comment='#', header=None,
                         usecols=(0,3,5,15,16), names=['orf', 'gene', 'qlen', 'qstart', 'qend'])
    if not len(data):
        return []
    data['gene'] = data['gene'].map(lambda m: normalize_marker_trans__dict.get(m , m))
    qlen = data[['gene','qlen']].drop_duplicates().set_index('gene')['qlen']

    def contig_name(ell):
        if orf_finder in ['prodigal', 'fast-naive']:
            contig,_ = ell.rsplit( '_', 1)
        else:
            contig,_,_,_ = ell.rsplit( '_', 3)
        return contig

    data = data.query('(qend - qstart) / qlen > 0.4').copy()
    data['contig'] = data['orf'].map(contig_name)
    contig_set = set(contig_names)
    data = data[data.apply(lambda row: row['contig'] in contig_set, axis=1)]
    if min_contig_len is not None:
        contig_len = {h:len(seq) for h,seq in fasta_iter(fasta_path)}
        data = data[data['contig'].map(lambda c: contig_len[c] >= min_contig_len)]
    data = data.drop_duplicates(['gene', 'contig'])

    from collections import defaultdict
    marker = data['gene'].values
    contig = data['contig'].values
    sequence2markers = defaultdict(list)
    for m, c in zip(marker, contig):
        sequence2markers[c].append(m)
    return sequence2markers


def generate_markers(fasta_path, binned_length, num_process, output = None, orf_finder = 'prodigal', prodigal_output_faa=None):
    '''Estimate number of bins from a FASTA file

    Parameters
    fasta_path: path
    binned_length: int (minimal contig length)
    num_process: int (number of CPUs to use)
    '''
    from .orffinding import run_orffinder
    with tempfile.TemporaryDirectory() as tdir:
        if output is not None:
            if os.path.exists(os.path.join(output, 'markers.hmmout')):
                return get_marker(os.path.join(output, 'markers.hmmout'), fasta_path, binned_length, orf_finder=orf_finder)
            else:
                os.makedirs(output, exist_ok=True)
                target_dir = output
        else:
            target_dir = tdir

        contig_output = run_orffinder(fasta_path, num_process, tdir, orf_finder, prodigal_output_faa=prodigal_output_faa)

        hmm_output = os.path.join(target_dir, 'markers.hmmout')
        try:
            with open(os.path.join(tdir, 'markers.hmmout.out'), 'w') as hmm_out_log:
                subprocess.check_call(
                    ['hmmsearch',
                     '--domtblout',
                     hmm_output,
                     '--cut_tc',
                     '--cpu', str(num_process),
                     os.path.split(os.path.abspath(__file__))[0] + '/marker.hmm',
                     contig_output,
                     ],
                    stdout=hmm_out_log,
                )
        except:
            if os.path.exists(hmm_output):
                os.remove(hmm_output)
            sys.stderr.write(
                f"Error: Running hmmsearch fail\n")
            sys.exit(1)

        return get_marker(hmm_output, fasta_path, binned_length, orf_finder=orf_finder)


def process_fasta(fasta_path):
    """
    Returns

    contigs: dictionary ID -> contig sequence
    """
    whole_contig_bp = 0
    contig_length_list = []
    contig_dict = {}

    for h, seq in fasta_iter(fasta_path):
        contig_length_list.append(len(seq))
        whole_contig_bp += len(seq)
        contig_dict[h] = seq

    return contig_dict

