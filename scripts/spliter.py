from Bio import SeqIO
#import argparse

def spliter(msa, keys, idx):
    msa_as_list = [[x.id, str(x.seq)] for x in SeqIO.parse(msa, 'fasta')]
    keysidx = dict(zip(idx, keys))

    for entry in msa_as_list:   
        idx = keysidx[entry[0]]
        with open(f'./aligned/{idx}.fasta', 'w') as fasta:
	        fasta.write(f'>{idx}\n{entry[1]}\n')

if __name__ == '__main__':
    msa = str(snakemake.input)
    keys = list(snakemake.params.k)
    ids = list(snakemake.params.i)
    
    spliter(msa, keys, ids)
