import numpy as np
from Bio import SeqIO
import argparse
import subprocess
import pandas as pd

def splits_for_blast(target, NAME):
    target_seq = target[1:]
    no_gap_t = target_seq[target_seq != '-']
    target_length = no_gap_t.shape[0]
    
    if target_length < 400:
        return 'less than 400 nucleotides'
    else:
        sub_aligns =[ [[f'>{NAME}_{x}'] , no_gap_t[x:x+400]] for x in range(0, target_length, 50) if  len(no_gap_t[x:x+400]) == 400]
        return sub_aligns


def sub_aligns(target, NAME):
    seq_target = [np.array([j for j in [x.id]] + list(str(x.seq))) for x in SeqIO.parse(target, 'fasta')][0]  
    sub_alignments = splits_for_blast(seq_target, NAME)


    if sub_alignments == 'less than 400 nucleotides':
        return 'impossible to test recomb (lenght < 400bp)'
    else:
        name_out = f'blast/sub_{NAME}.fasta'
        
        with open(name_out, 'w') as out_subs:
            for sub in sub_alignments:
                out_subs.write(f'{sub[0][0]}\n{"".join(sub[1])}\n')

        return name_out


def do_blast_window(target, NAME):
    to_blast = sub_aligns(target, NAME)
    if to_blast == 'impossible to test recomb (lenght < 400bp)':
        with open(f'blast/recblast_{NAME}.txt', 'w') as out_subs:
            out_subs.write('impossible to test recomb (lenght < 400bp)\n')
    else:
        subprocess.call(['blastn', '-db', './data/db_01-02_and_pures', '-query',
        '{0}'.format(to_blast), '-out', 'blast/recblast_{0}.txt'.format(NAME),
        '-word_size', '10', '-outfmt', '10', '-evalue', '1.e-100'])

        df = pd.read_csv(f'blast/recblast_{NAME}.txt', header=None)
        df[[0,1,10,11]].to_csv(f'blast/recblast_{NAME}.txt', index=False)        

    comand_rm = 'rm {0}'.format(to_blast)
    subprocess.call(comand_rm.split(' '))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", required = True,
                            help="Name input file.")
    args = parser.parse_args()
    NAME = str(args.input).replace('aligned/aligned_' ,'').replace('.fasta', '')
    do_blast_window(args.input, NAME)
