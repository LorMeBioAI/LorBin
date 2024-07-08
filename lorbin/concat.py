from Bio import SeqIO
import argparse
import pandas as pd

if __name__=="__main__":
    parser= argparse.ArgumentParser(description="""Creates the input FASTA file for LorBin.
    Input should be one or more FASTA files, each from a sample-specific assembly.
    resulting FASTA can be binsplit with separator '-'.""",
    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o","--output",help="Path to output FASTA file", required=True)
    parser.add_argument("-fa","--fasta",help="Paths to input FASTA files", nargs="+", required=True)

    args = parser.parse_args()
    fastas = args.fasta
    fw = open(f'{args.output}','w')
    for e,fasta in enumerate(fastas):
        for record in SeqIO.parse(fasta, 'fasta'):
            fw.write(f'>{e}-{record.id}\n')
            fw.write(f'{record.seq}\n')
    fw.close()
    #pd.DataFrame(fastas).to_csv(f"{args.output.split('.')[0]}.csv")

