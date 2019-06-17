.. _usage:

SNAPPy commands
===============

This is a list of the implemented functionalities in SNAPPy:

- Basic subtyping::

    snakemake subtype_all

- Subtyping using n cpu cores::

    snakemake subtype_all --cores n


- Basic alignment::

    snakemake align_all


- Alignment using n cpu cores::

    snakemake align_all --cores n


- Clean the SNAPPy folder from previous runs outputs::

    snakemake delete_all_outputs

- Clean the SNAPPy folder from previous runs intermediate files (but not outputs)::

    snakemake delete_interm_files

- Compress to a file and clean the SNAPPy folder from previous runs intermediate files (but not outputs)::

    snakemake compress_and_del_inter_files

- Compress everything inside the SNAPPy folder to a file::

    snakemake compress_all

Note: Please note that SNAPPy is built on top of `Snakemake <https://snakemake.readthedocs.io/en/stable/index.html>`_, and therefore should retain all the functionalities of this pipeline software.
