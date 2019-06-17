import argparse
import subprocess
from Bio import SeqIO
import pandas as pd
import numpy as np

ALL_REFS = [[x.id] + list(str(x.seq)) for x in SeqIO.parse(f'data/all_refs.fasta', 'fasta')]
OUT_CPZ = [[x.id] + list(str(x.seq)) for x in SeqIO.parse(f'data/out_cpz.fasta', 'fasta')][0]

def blast_closser(file_name, NAME):

    out_nogap = open('blast/nogap_{0}.fasta'.format(NAME), "w")
    subprocess.call(['sed', '-e', 's/-//g', '{0}'.format(file_name)], stdout=out_nogap)

    subprocess.call(['blastn', '-db', 'data/db_all_refs', '-query',
    'blast/nogap_{0}.fasta'.format(NAME), '-out', 'blast/blast_{0}.txt'.format(NAME),
    '-word_size', '10', '-outfmt', '10', '-evalue', '1.e-10'])

    subprocess.call(['rm', 'blast/nogap_{0}.fasta'.format(NAME)])

    try:
        df = pd.read_csv(f'blast/blast_{NAME}.txt', header=None)    
        filter_df = df.sort_values(by=[11], ascending=False).copy()
        
        top_48_closser = filter_df[1].values[:48]
        top_48_pures = [x for x in filter_df[1].values if ('_' not in x) | (x[:2] in ['01', '02'])][:48]
        top_48_recombs = [x for x in filter_df[1].values if ('_' in x) & (x[:2] not in ['01', '02'])][:48]
        
        df[[0,1,10,11]].to_csv(f'blast/blast_{NAME}.txt', index=False)

        return [NAME, top_48_closser, top_48_pures, top_48_recombs]
        
    except:
        print(f'ERROR while performing BLAST of {NAME}')
        return [NAME, [], [], []]
        

def build_msas(file_name, NAME):
    closser_data = blast_closser(file_name, NAME)
    target = [np.array([x.id] + list(str(x.seq))) for x in SeqIO.parse(f'{file_name}', 'fasta')][0]
    mask_gaps = target != '-'

    specify_type_msas(1, 'all', target, mask_gaps, closser_data, NAME)
    specify_type_msas(2, 'pure', target, mask_gaps, closser_data, NAME)
    specify_type_msas(3, 'recomb', target, mask_gaps, closser_data, NAME)


def specify_type_msas(col, name, target, mask_gaps, closser_data, NAME):
    msa_idx = closser_data[col]
    closser_seqs = np.array([x for x in ALL_REFS if x[0] in msa_idx])

    masked_target = target[mask_gaps]
    masked_root = np.array(OUT_CPZ)[mask_gaps].copy()
    masked_refs = closser_seqs[:, mask_gaps].copy()
 
    with open(f'trees/msa_{name}_{NAME}.fasta', 'w') as out_fasta:
        for seq in masked_refs:
            out_fasta.write(f'>{seq[0]}\n{"".join(seq[1:])}\n')
                            
        out_fasta.write(f'>{masked_target[0]}\n{"".join(masked_target[1:])}\n')
        out_fasta.write(f'>{masked_root[0]}\n{"".join(masked_root[1:])}\n')

def tree_maker(file_name, NAME):
    build_msas(file_name, NAME)

    out_all = open("trees/all_{0}.nwk".format(NAME), "w")
    subprocess.call(['fasttree', '-quiet', '-gtr', '-nopr', '-nt',
    'trees/msa_all_{0}.fasta'.format(NAME)], stdout=out_all, stderr=subprocess.STDOUT)
    subprocess.call(['rm', 'trees/msa_all_{0}.fasta'.format(NAME)])

    out_pure = open("trees/pure_{0}.nwk".format(NAME), "w")
    subprocess.call(['fasttree', '-quiet', '-gtr', '-nopr', '-nt',
    'trees/msa_pure_{0}.fasta'.format(NAME)], stdout=out_pure, stderr=subprocess.STDOUT)
    subprocess.call(['rm', 'trees/msa_pure_{0}.fasta'.format(NAME)])

    out_recomb = open("trees/recomb_{0}.nwk".format(NAME), "w")
    subprocess.call(['fasttree', '-quiet', '-gtr', '-nopr', '-nt',
    'trees/msa_recomb_{0}.fasta'.format(NAME)], stdout=out_recomb, stderr=subprocess.STDOUT)
    subprocess.call(['rm', 'trees/msa_recomb_{0}.fasta'.format(NAME)])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", required = True,
                            help="Name input file.")
    parser.add_argument("-g","--genomic", required = True,
                            help="Name of the genomic region.")

    args = parser.parse_args()
    NAME = str(args.input).replace('aligned/aligned_' ,'').replace('.fasta', '')
    
    gr = args.genomic
    if gr != 'GAG-POL-ENV':
        print('GENOMIC REGION IN CONFIG FILE IS WRONG! FOR SUBTYPING USE ONLY :"GAG-POL-ENV"\n')
        exit()
    else:
        pass

    tree_maker(args.input, NAME)
