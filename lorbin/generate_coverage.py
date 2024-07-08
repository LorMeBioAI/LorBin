import os
import subprocess
from .atomicwrite import atomic_write

def calculate_coverage(depth_stream, bam_file, edge=75,contig_threshold=1000):
    """
    depth_stream : an iterable like the output of bedtools genomecov
    bam_file : str filename
    edge : int

    """
    import pandas as pd
    import numpy as np
    from itertools import groupby

    contigs = []
    mean_coverage = []

    for contig_name, lines in groupby(depth_stream, lambda ell: ell.split('\t', 1)[0]):
        lengths = []
        values = []
        for line in lines:
            line_split = line.strip().split('\t')
            length = int(line_split[2]) - int(line_split[1])
            value = int(float(line_split[3]))
            lengths.append(length)
            values.append(value)
        depth_value = np.zeros(sum(lengths), dtype=int)
        s = 0
        for ell,v in zip(lengths, values):
            depth_value[s:s+ell] = v
            s += ell

        cov_threshold = contig_threshold
        if len(depth_value) < cov_threshold:
            continue
        depth_value_ = depth_value[edge:-edge]
        mean_coverage.append(depth_value_.mean())
        contigs.append(contig_name)


    contig_cov = pd.DataFrame(
            { '{0}_cov'.format(bam_file): mean_coverage,
             }, index=contigs)
    return contig_cov


def generate_cov(bam_file, bam_index, out, logger, contig_threshold=1000):
    """
    Call bedtools and generate coverage file

    bam_file: bam files used
    out: output
    threshold: threshold of contigs that will be binned
    sep: separator for multi-sample binning
    """
    import numpy as np
    logger.debug('Processing `{}`'.format(bam_file))
    bam_name = os.path.split(bam_file)[-1] + '_{}'.format(bam_index)

    bed_p = subprocess.Popen(
        ['bedtools', 'genomecov',
         '-bga',
         '-ibam', bam_file],
        universal_newlines=True,
        stdout=subprocess.PIPE)

    contig_cov = calculate_coverage(bed_p.stdout, bam_file, contig_threshold=contig_threshold)
    if bed_p.wait() != 0:
        logger.critical(f"Running `bedtools genomecov` failed ({bam_file}). Please check your input files: LorBin expects that they are sorted BAM files.")
        raise OSError(f"Failure running `bedtools genomecov` ({bam_file})")
    elif len(contig_cov) == 0:
        logger.critical("Running `bedtools genomecov` did not return an error, but the result is an empty file. Please check your input files: LorBin expects that they are sorted BAM files.")
        raise OSError("Running bedtools returned an empty file")

    contig_cov = contig_cov.apply(lambda x: x + 1e-5)
    with atomic_write(os.path.join(out, '{}_data_cov.csv'.format(bam_name)), overwrite=True) as ofile:
        contig_cov.to_csv(ofile)

    return bam_file

def combine_cov(cov_dir, bam_list):
    """
    generate cov for every sample in one file
    """
    import pandas as pd
    data_cov = pd.read_csv(os.path.join(cov_dir, '{}_data_cov.csv'.format(
        os.path.split(bam_list[0])[-1] + '_{}'.format(0))), index_col=0)

    for bam_index, bam_file in enumerate(bam_list):
        if bam_index == 0:
            continue
        cov = pd.read_csv(os.path.join(cov_dir, '{}_data_cov.csv'.format(
                os.path.split(bam_file)[-1] + '_{}'.format(bam_index))),index_col=0)
        cov.index = cov.index.astype(str)
        data_cov = pd.merge(data_cov, cov, how='inner', on=None,
                            left_index=True, right_index=True, sort=False, copy=True)

    return data_cov

