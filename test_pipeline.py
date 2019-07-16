from scripts import cleaner
from scripts import spliter
from scripts import align_and_map
from scripts import concatenator
from scripts.tree_maker import *
from scripts.blast_recomb import *
from scripts.summarize_results import *

import os
import subprocess


KEYS = ['1', '2', '3']
IDX = ['test01', 'test02', 'test03']
GR = 'GAG-POL-ENV'

# Test config
def test_config():
    with open('config.yaml', 'r') as read_config:
        config = read_config.readlines()

    assert config[1].split("'")[1] == GR

# Test cleaner
## Create test folder
directory = 'test/testfolder'
if not os.path.exists(directory):
    os.makedirs(directory)

## Create test file
open('test/testfile', 'a').close()

## Tests
def test_clean_folders():
    cleaner.clean_folders('test/testfolder')
    assert os.path.exists('test/testfolder') == False

def test_clean_files():
    cleaner.clean_files('test/testfile')
    assert os.path.exists('test/testfile') == False

# Test spliter
## Create test folder
directory = 'aligned'
if not os.path.exists(directory):
    os.makedirs(directory)

## Tests
def test_spliter():
    spliter.spliter('test/test_msa.fasta', KEYS, IDX)
    assert sorted(os.listdir('aligned')) == sorted(['1.fasta', '2.fasta', '3.fasta'])

# Test align and map
## Tests
def test_align_and_map():
    inputs = os.listdir('aligned')
    expected_output = ['aligned_1.fasta', 'aligned_2.fasta', 'aligned_3.fasta']

    for file in inputs:
        align_and_map.align_and_map(f'aligned/{file}', f'aligned/aligned_{file}', GR)

    real_output = [x for x in os.listdir('aligned') if 'aligned' in x]
    assert sorted(real_output) == sorted(expected_output)

# Test concatenator
## Tests
def test_concatenator():
    inputs = [x for x in sorted(os.listdir('aligned')) if 'aligned' in x]

    with open('all_aligned.fasta', 'w') as out_msa:
        for file in inputs:
            seq = concatenator.concatenator(f'aligned/{file}', KEYS, IDX)
            out_msa.write(f'>{seq[0]}\n{seq[1]}\n')

    with open('all_aligned.fasta') as this_output:
        out = this_output.read()
    
    with open('test/all_aligned.fasta') as expected_output:
        e_out = expected_output.read() 

    assert out == e_out


# Test tree_maker
## Create test folders
directory = 'blast'
if not os.path.exists(directory):
    os.makedirs(directory)

directory = 'trees'
if not os.path.exists(directory):
    os.makedirs(directory)

## Tests blast_closser
def test_blast_closse():
    inputs = [x for x in os.listdir('aligned') if 'aligned' in x]
    
    for file in inputs:
        NAME = str(f'aligned/{file}').replace('aligned/aligned_' ,'').replace('.fasta', '')
        blast_closser(f'aligned/{file}', NAME)

    expected_output = ['blast_1.txt', 'blast_2.txt', 'blast_3.txt']
    real_output = [x for x in os.listdir('blast') if 'blast_' in x]
    assert sorted(real_output) == sorted(expected_output)

## Tests build_msas
def test_build_msas():
    inputs = [x for x in os.listdir('aligned') if 'aligned' in x]

    for file in inputs:
        NAME = str(f'aligned/{file}').replace('aligned/aligned_' ,'').replace('.fasta', '')
        build_msas(f'aligned/{file}', NAME)

    expected_output = ['msa_recomb_2.fasta', 'msa_all_2.fasta', 'msa_pure_3.fasta',
    'msa_all_3.fasta', 'msa_recomb_3.fasta', 'msa_all_1.fasta', 'msa_pure_1.fasta',
    'msa_pure_2.fasta', 'msa_recomb_1.fasta']

    real_output = [x for x in os.listdir('trees') if 'msa_' in x]

    assert sorted(real_output) == sorted(expected_output)

## Tests tree_maker
def test_tree_maker():
    inputs = [x for x in os.listdir('aligned') if 'aligned' in x]

    for file in inputs:
        NAME = str(f'aligned/{file}').replace('aligned/aligned_' ,'').replace('.fasta', '')
        tree_maker(f'aligned/{file}', NAME)

    expected_output = ['pure_2.nwk', 'pure_1.nwk', 'all_2.nwk', 'recomb_3.nwk',
    'all_1.nwk', 'recomb_2.nwk', 'all_3.nwk', 'pure_3.nwk', 'recomb_1.nwk']
    real_output = [x for x in os.listdir('trees') if '.nwk' in x]

    assert sorted(real_output) == sorted(expected_output)


# Test blast_recomb
## Tests do_sub_aligns
def test_do_sub_aligns():
    inputs = [x for x in os.listdir('aligned') if 'aligned' in x]
    for file in inputs:
        NAME = str(f'aligned/{file}').replace('aligned/aligned_' ,'').replace('.fasta', '')
        do_sub_aligns(f'aligned/{file}', NAME)

    expected_output = ['sub_2.fasta', 'sub_3.fasta', 'sub_1.fasta']
    real_output = [x for x in os.listdir('blast') if 'sub' in x]

    assert sorted(real_output) == sorted(expected_output)

## Tests do_blast_window
def test_do_blast_window():

    inputs = [x for x in os.listdir('aligned') if 'aligned' in x]
    for file in inputs:
        NAME = str(f'aligned/{file}').replace('aligned/aligned_' ,'').replace('.fasta', '')
        do_blast_window(f'aligned/{file}', NAME)

    expected_output = ['recblast_3.txt', 'recblast_1.txt', 'recblast_2.txt']
    real_output = [x for x in os.listdir('blast') if 'rec' in x]

    assert sorted(real_output) == sorted(expected_output)

# Test summarize_results
## Tests
def get_clossest_blast():

    all_trees_inputs = [x for x in os.listdir('trees') if (x[-4:]  == '.nwk') & (x[:4] == 'all_')]
    pure_trees_inputs = [x for x in os.listdir('trees') if (x[-4:]  == '.nwk') & (x[:5] == 'pure_')]
    recomb_trees_inputs = [x for x in os.listdir('trees') if (x[-4:]  == '.nwk') & (x[:7] == 'recomb_')]
    blast_c = [x for x in os.listdir('blast') if (x[-4:]  == '.txt') & (x[:6] == 'blast_')]
    blast_inputs = [x for x in os.listdir('blast') if (x[-4:]  == '.txt') & (x[:9] == 'recblast_')]

    results = {}
    for pos, key in enumerate(KEYS):
        results[key] = [IDX[pos]]

    for blast_res in blast_inputs:
        output = process_blast_recomb(blast_res)
        results[output[0]] += output[1]

    for tree in all_trees_inputs:
        output = process_trees(tree)
        results[output[0]] += output[1:]
    
    for tree in pure_trees_inputs:
        output = process_trees(tree)
        results[output[0]] += output[1:]

    with open('subtype_results.csv') as this_output:
        out = this_output.read()
    
    with open('test/subtype_results.csv') as expected_output:
        e_out = expected_output.read() 

    with open('report_subtype_results.csv') as this_output:
        rout = this_output.read()
    
    with open('test/report_subtype_results.csv') as expected_output:
        e_rout = expected_output.read() 


    assert (out == e_out) & (rout == e_rout)

# Test cleaner all
## Tests
def test_delete_all_outputs():
    FOLDERS = ['trees','aligned','blast']
    FILES = ['keys_and_ids.csv','subtype_results.csv', 'report_subtype_results.csv.csv',
    'all_aligned.fasta']
    for folder in FOLDERS:
        cleaner.clean_folders(folder)
    for file in FILES:
        cleaner.clean_files(file)


