from Bio import SeqIO

def spliter(msa, keys, idx):
    """Spliter

    This function uses the Biopython module SeqIO to parse a multiple sequences 
    alignment and separate each fasta sequence in a separated file. The
    resulting files are written to the folder 'aligned' with the following
    notation: {id_of_the_fasta_sequence}.fasta.

	Args:
        msa (fasta): Multiple sequences alignment in fasta format.
        keys (str): List of fasta headers in the msa file.
        idx (str): List of internal SNAPPy ids corresponding to the header in 
        the msa file.

	Returns:
        This function does not return.
	"""
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
