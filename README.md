# Snakemake subtyper

This is snakemake pipeline represents my attemp of automatizing the alignment of HIV-1 fasta sequences, for the purpose of phylogenetic analyses, with the the possibility of performing a subtype analysis based on phylogenetic pairing with subtype reference sequences.

This was designed to execute in unix based systems. However, as long as one can install all the dependencies it should work. To facilitate the instalation of dependencies i strongly advice the usage of anaconda or miniconda as a package manager.

## How to install

1. Strongly advise the usage of conda (anaconda or miniconda) as a package manager (instalation instructions here: https://conda.io/docs/user-guide/install/index.html). The following steps assume a proper instalation of conda in the system;
2. Create the environment, this should one be done once per system: conda-env create -f environment.yaml;
3. Confirm if the environment was created/you already had it: conda-env list;

    3.1. If "subtype_hiv1" is in the list good, otherwise go to step 2;
4. Activate the environmet:  source activate subtype_hiv1 ( works in unix systems);
5. The snakemake pipline is ready to use, please see the next topic on the running options.

## How to use

1. Place your multiple sequence file in fasta format in the input folder;
2. Open the "config.yaml" file and write the genomic region you want  to align into the line under "genomic region:" (Valid inputs:'PR', 'RT', 'PR-RT', 'INT');
3. To perform the multiple sequence aligment run: "snakemake align_all"

    3.1. Multiple sequence alignment will be outputed as : "all_aligned.fasta"
4. To perform subtyping run : "snakemake subtype_all";

    4.1. Please note that for subtyping the sequences need to be previously aligned, the pipeline will understand this need and perform it if needed;

    4.2. Two resulting files will be outputed; "all_subtyped.csv": a CSV file (separated by ';') with the sequence name and the inferred subtype; and "report_all_subtyped.csv": a csv file (separated by ';') with the sequence name and information about the closses nodes used to perform the subtyping.