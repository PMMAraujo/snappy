import argparse
import subprocess
from Bio import SeqIO
import pandas as pd
import numpy as np

ALL_REFS = [[x.id] + list(str(x.seq)) for x in SeqIO.parse(f'data/all_refs.fasta', 'fasta')]
OUT_CPZ = [[x.id] + list(str(x.seq)) for x in SeqIO.parse(f'data/out_cpz.fasta', 'fasta')][0]

def blast_closser(file_name, NAME):
    """BLAST to find closser reference

    From a single fasta file a version of it without any gaps ('-') is created
    and written to the folder 'blast' as 
    'nogap_{id_of_the_fasta_sequence}.fasta'.
    Then that no gap sequence is used as an input in a BLAST against the
    database 'data/db_all_refs'.
    Three group of reference sequences are created, those will be used fot three
    diferent phylogenetic analysis. One contains the first 48 BLAST results,
    the 'pures' contains the first 48 results excluding circulating recombinat
    forms (CRFs), and the 'recomb' contains the first 48 results for only CRFs.

    This function calls MAFFT to performe a sequence alignment between the
    HIV-1 reference sequence (NCBI id: HXB2) and a sequence in a given file
    'msa'. The alignment method does not allow any gaps the reference sequence.
    After the alignment is performed the given sequence in 'msa' will be trimed
    to only contain the genomic region specified by the user in 'gr'. The
    resulting file is them written to the folder 'aligned' with the following
    notation: aligned_{id_of_the_fasta_sequence}.fasta. This function is called
    by the function 'build_msas'.

	Args:
        file_name (fasta): Fasta file containing 1 fasta sequence aligned.
        NAME (str): Global variable. Internal index of SNAPPy for this fasta.

	Returns:
        List of lists of reference sequences names.
	"""
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
    """Creating the MSA for phylogenetic inference

    This is just a wraper function on top of the function 'specify_type_msas'
    to execute it three times.
    This function calls the functions 'blast_closser' and 'specify_type_msas' 
    and is called by the function 'tree_maker'.

	Args:
        file_name (fasta): Fasta file containing 1 fasta sequence aligned.
        NAME (str): Global variable. Internal index of SNAPPy for this fasta.

	Returns:
        This function does not return.
	"""
    closser_data = blast_closser(file_name, NAME)
    target = [np.array([x.id] + list(str(x.seq))) for x in SeqIO.parse(f'{file_name}', 'fasta')][0]
    mask_gaps = target != '-'

    specify_type_msas(1, 'all', target, mask_gaps, closser_data, NAME)
    specify_type_msas(2, 'pure', target, mask_gaps, closser_data, NAME)
    specify_type_msas(3, 'recomb', target, mask_gaps, closser_data, NAME)


def specify_type_msas(col, name, target, mask_gaps, closser_data, NAME):
    """Creating the MSA for phylogenetic inference


    This writtes msa files. From the BLAST output it infers the best result
    (bitscore), geting all references with that score. Then the positions with
    gaps present in the target sequence are deleted in the outgroup sequence
    ('OUT_CPZ') and the remaining closser references. Finaly the msa file is
    written to the folder 'trees' with the following notation:
    msa_{type}_{id_of_the_fasta_sequence}.fasta. Type can take the any string
    but only 'all, 'pure', and 'recomb' are used.

	Args:
        col (int): Number of the column where the list of selected referes is.
        name (str): Can take the any string but only 'all, 'pure', and 'recomb'
        are used.
        target (np.array): Target fasta in an array.
        mask_gaps (Bolean): Locations in the target sequence with gaps.
        closser_data (list): Output from the blast_closser function.
        NAME (str): Global variable. Internal index of SNAPPy for this fasta.

	Returns:
        This function does not return.
	"""
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
    """Phylogenetic inference


    This function uses the multiple sequence alignments created by the function 
    'build_msas' using the function 'specify_type_msas'.
    Three trees are created per sequence 'id' (file_name). The trees are
    outputed to the folder 'trees' with the following notation:
    {type}_{id_of_the_fasta_sequence}.nwk. Type can take the any string
    but only 'all, 'pure', and 'recomb' are used.

	Args:
        file_name (str): Name of the input fasta file.
        NAME (str): Global variable. Internal index of SNAPPy for this fasta.

	Returns:
        This function does not return.
	"""
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

    is_empty = str(list(SeqIO.parse(args.input, 'fasta'))[0].seq).replace('-', '')
    if is_empty == '': # deal with empty fastas after alignment
        with open('blast/blast_{0}.txt'.format(NAME), "w") as out_b:
            out_b.write('not enough genomic information\n')
        with open("trees/all_{0}.nwk".format(NAME), "w") as out_t:
            out_t.write('not enough genomic information\n')
        with open("trees/pure_{0}.nwk".format(NAME), "w") as out_t:
            out_t.write('not enough genomic information\n')
        with open("trees/recomb_{0}.nwk".format(NAME), "w") as out_t:
            out_t.write('not enough genomic information\n')
    else: # process sequences normaly
        tree_maker(args.input, NAME)
