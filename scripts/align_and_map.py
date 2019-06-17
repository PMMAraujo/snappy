import argparse
import subprocess


REGIONS = {'GAG':[789, 2292], 'PR':[2252, 2550],
		   'RT':[2549,3870], 'PR-RT':[2252,3869],
		   'INT':[4229,5096], 'POL':[2252, 5096],
		   'ENV':[6224, 8795],
		   'GAG-POL-ENV':[789, 2292, 2252, 5096, 6224, 8795] }

def align_and_map(msa, out, gr):
    """Align and Map.

    This function calls MAFFT to performe a sequence alignment between the
    HIV-1 reference sequence (NCBI id: HXB2) and a sequence in a given file
    'msa'. The alignment method does not allow any gaps the reference sequence.
    After the alignment is performed the given sequence in 'msa' will be trimed
    to only contain the genomic region specified by the user in 'gr'. The
    resulting filke is them wirtte to the folder 'aligned' with the following
    notation: aligned_{name_of_the_fasta_sequence}.fasta.

	Args:
        msa (fasta): Fasta file containing 1 fasta sequence to align.
        out (fasta): Name of the aligned version of the fasta file.
        gr (str): String of the genomic region of interest. Specified in the
        'config.yaml' file.

	Returns:
        This function doe not return.
	"""
    try:
        this_region = REGIONS[gr]
    except KeyError:
        print(f"Invalid Genomic regioin in the config file!!!\n(Valid \
        inputs: {list(REGIONS.keys())})\n")

    comand = 'mafft --quiet --add {0} --keeplength ./data/HXB2.fasta'.format(msa)
    pre = subprocess.check_output(comand.split(' '), universal_newlines=True)
    as_list = pre.split('\n')[:-1]
    id_place = [x for x in enumerate(as_list) if x[1][0] == '>'][1][0]
    seq_id = as_list[163][1:]
    all_seq = ''.join(as_list[int(id_place) +1:])

    if len(this_region) == 2:
        sequence = all_seq[this_region[0]:this_region[1]]
    elif len(this_region) == 6:
        gag_reg = all_seq[this_region[0]:this_region[1]]
        pol_reg = all_seq[this_region[2]:this_region[3]]
        env_reg = all_seq[this_region[4]:this_region[5]]
        sequence = gag_reg + pol_reg + env_reg
    else:
        print('ERROR while passing genomic region. Please check the \
        integrety of the script: align_and_map.')

    with open(out, 'w') as mapped:
        mapped.write('>{}\n{}\n'.format(seq_id, sequence))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", required = True,
                            help="Name input file.")
    parser.add_argument("-o","--output", required = True,
                            help="Name output file.")
    parser.add_argument("-g","--genomic", required = True,
                            help="Name of the genomic region.")
    args = parser.parse_args()  
    align_and_map(args.input, args.output, args.genomic)
