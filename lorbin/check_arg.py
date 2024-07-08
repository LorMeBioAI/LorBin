import os
import logging
def check_generate_data(logger, fastadir, bams, outdir):
    if fastadir.split('.')[-1] not in ['fa','fasta','fna']:
        logger.info('the format of input file must be .fa, .fasta, .fna')
        return False
    elif not os.path.exists(fastadir):
        logger.info(f'the {fastadir} not exists')
        return False
    for bam in bams:
        logger.info(bam)
        if bam.split('.')[-1] !='bam':
            logger.info('the format of --bams must be .bam')
            return False
        elif not os.path.exists(bam):
            logger.info(f'the {bam} not exists')
            return False
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        logger.info(f'the {outdir} makes')
    return True

def check_cluster(logger, outdir, fastadir, embeddingdir, feature, a):
    if fastadir.split('.')[-1] not in ['fa','fasta','fna']:
        logger.info('the format of input file must be .fa, .fasta, .fna')
        return False
    elif not os.path.exists(fastadir):
        logger.info(f'the {fastadir} not exists')
        return False
    if not os.path.exists(embeddingdir):
        return False
    if feature not in ['no_markers', 'markers110', 'markers35']:
        logger.info("not find th evaluation model, feature should be 'no_markers', 'markers110', 'markers35'.")
        return False
    if a>1 or a<0:
        logger.info('parameter -a must between 0 and 1')
        return False
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        logger.info(f'the {outdir} makes')
    return True

def check_train(logger, outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        logger.info(f'the {outdir} makes')
    return True

if __name__=="__main__":
    logger = logging.getLogger('LorBin')
    print(check_generate_data(logger,'../test/CRR451057.hifiasm.fna',['../test/CRR451057.sorted.mapped.bam'],'../test_o'))
