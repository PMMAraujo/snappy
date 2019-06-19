### IMPORTS AND DEPENDENCIES
import os
from Bio import SeqIO
import sys
import subprocess
from scripts import keys_creator
import shutil


configfile: "config.yaml"


### GLOBAL VARIABLES

# define genomic region from config file
GR = config['genomic_region']

# confirm there is only 1 file in input folder and define it as INPUT
if len(os.listdir('input')) > 1:
    sys.exit('ERROR: more than 1 file in the input folder')
else:
    INPUT = os.listdir('input')[0]

# get ids from the INPUT file (which is as msa)
MSA = SeqIO.parse('input/{}'.format(INPUT), 'fasta')
KEYS, IDS = keys_creator.create_keys([seq.id for seq in MSA])


### SUBTYPER

rule subtype_all:
    input:
        trees_all=expand("trees/all_{sample}.nwk", sample=KEYS),
        trees_pure=expand("trees/pure_{sample}.nwk", sample=KEYS),
        trees_recomb=expand("trees/recomb_{sample}.nwk", sample=KEYS),
        blast_c=expand("blast/blast_{sample}.txt", sample=KEYS),
        blast=expand("blast/recblast_{sample}.txt", sample=KEYS)
    output:
        "subtype_results.csv"
    params:
        k=lambda wildcards: [x for x in KEYS],
        i=lambda wildcards: [x for x in IDS]
    script:
        "scripts/summarize_results.py"

### BLAST RECOMBINATION
rule blast_recomb:
    input:
        "aligned/aligned_{sample}.fasta"
    output:
        "blast/recblast_{sample}.txt"
    shell:
        "python scripts/blast_recomb.py -i '{input}'"


### TREE BUILDING

rule tree_maker:
    input:
        "aligned/aligned_{sample}.fasta"
    output:
        "trees/all_{sample}.nwk",
        "trees/pure_{sample}.nwk",
        "trees/recomb_{sample}.nwk",
        "blast/blast_{sample}.txt"
    priority:
        50
    shell:
        "python scripts/tree_maker.py -i '{input}' -g '{GR}'"


### ALIGNER

rule align_all:
    input:
        files = expand("aligned/aligned_{sample}.fasta", sample=KEYS)
    output:
        "all_aligned.fasta"
    params:
        k=lambda wildcards: [x for x in KEYS],
        i=lambda wildcards: [x for x in IDS]
    script:
        "scripts/concatenator.py"


rule map_and_align:
    input:
        "aligned/{sample}.fasta"
    output:
        "aligned/aligned_{sample}.fasta"
    shell:
        "python scripts/align_and_map.py -i '{input}' -o '{output}' -g '{GR}'"


rule split:
    input:
        "input/{}".format(INPUT)
    output:
        expand("aligned/{sample}.fasta", sample=KEYS)
    params:
        k=lambda wildcards: [x for x in KEYS],
        i=lambda wildcards: [x for x in IDS]
    script:
        "scripts/spliter.py"


### HELPERS(quality of life)
onsuccess:
    shutil.rmtree(".snakemake")

rule delete_interm_files:
    params:
        "inter"
    shell:
        "python scripts/cleaner.py -w {params}"

rule delete_all_outputs:
    params:
        "all"
    shell:
        "python scripts/cleaner.py -w {params}"

rule compress_and_del_inter_files:
    params:
        "c_inter"
    shell:
        "python scripts/cleaner.py -w {params}"

rule compress_all:
    params:
        "c_all"
    shell:
        "python scripts/cleaner.py -w {params}"
