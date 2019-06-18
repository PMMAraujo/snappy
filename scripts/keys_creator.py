
def create_keys(ids):
    """Keys creator

    This function creates a list of integers to be used as internal SNAPPy ids
    to avioid issues with compex fasta file names.

	Args:
        ids (list): List of fasta headers in the msa file.

	Returns:
        List of pairs (keys, ids)
	"""
    d = len(str(len(ids)))
    new_keys = [f"%0{d}d" % x for x in  list(range(len(ids)))]

    return [new_keys, ids]
