from Bio import SeqIO

def concatenator(input, keys, ids):
    """Fasta files concatenator

    This function uses the Biopython module SeqIO to parse several single 
    sequences fata files. It then procedes to write a single multiple sequence
    alignment ( of aligned sequences).

	Args:
        input (fasta): Single sequence fasta file.
        keys (str): List of fasta headers in the msa file.
        idx (str): List of internal SNAPPy ids corresponding to the header in 
        the msa file.
	Returns:
        all_aligned.fasta ( fasta): Multiple sequence alignment file
	"""
    this_seq = [[x.id, str(x.seq)] for x in SeqIO.parse(input, 'fasta')][0]
    where = keys.index(this_seq[0])

    return [ids[where], this_seq[1]]

if __name__ == '__main__':
	files = list(snakemake.input)
	keys = list(snakemake.params.k)
	ids = list(snakemake.params.i)

	with open('all_aligned.fasta', 'w') as out_msa:

		for file in files:
			seq = concatenator(file, keys, ids)
			out_msa.write(f'>{seq[0]}\n{seq[1]}\n')
