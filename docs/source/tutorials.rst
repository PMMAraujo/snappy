.. _tutorials:

Tutorials
=========

Please ensure SNAPPy's installation was successful. If you need help regarding SNAPPy's installation please visit the :ref:`installation` section. To use SNAPPy the terminal always needs to be executed from the snappy folder, and the conda environment is active::

    conda actiavte snappy
 
If the tool is stopping abruptly or behaving unexpectedly ensure it is functional by simply running the tests::

    py.test

To run the tests the conda snappy environment need to be active. If the tests result in error it is advised to reinstall the tool.

 .. _tut1:

1: Subtyping with SNAPPy
^^^^^^^^^^^^^^^^^^^^^^^^

Disclaimer: This tutorial and Tutorial :ref:`tut2` share the first four steps.

HIV-1 subtyping is the main reason SNAPPY was build. To use SNAPPy the terminal always needs to be executed from the snappy folder, and the conda environment is active::

    conda actiavte snappy

We highly encourage to run the tests before using the pipeline just to make sure the tasks will run as intended and the pipeline was not modified or damaged::

    py.test

1)
Look at the content of the snappy folder::

    ls

Expected result::

    config.yaml  data  environment.yaml  input  LICENSE.txt  README.md  scripts  Snakefile  snappy_installer_x86_64.sh  test  test_pipeline.py

 
2)
In this tutorial we will focus on the folder 'input' and the  file 'config.yaml'. For exploration of the remaining content please look at the :ref:`how_it_works` section. Lets look inside the folder input::

    ls input

Expected result::

    test_msa.fasta

3)
As you can see SNAPPy has a file inside the folder 'input'. This is a test file used in this tutorial. For a successful run only one file can be inside the 'input' folder. Therefore, if you want to run an analysis with your data you need to first delete the default file inside the 'input' folder.

4)
Now lets look at the 'config.yaml' file:´

    cat config.yaml


Expected result::

    genomic_region:
      'GAG-POL-ENV'

5)
As you can see the genomic region selected for the SNAPPy's analysis is 'GAG-POL-ENV' which basically consists the HIV-1 regions named GAG, POL, and ENV concatenated (according to `HXB2 (K03455) <https://www.hiv.lanl.gov/components/sequence/HIV/asearch/query_one.comp?se_id=K03455>`_ reference genome). To use SNAPPy for subtyping the selected genomic region in the 'config.yaml' file always needs to be 'GAG-POL-ENV'. If for any reason you have changed it (for instance during Tutorial :ref:`tut2`) change it back.

6)
In order to the subtyping simple type::

    snakemake subtype_all

7)
The SNAPPy pipeline will then be executed and at green you should see each task being done. It should terminate with the following line::

    11 of 11 steps (100%) done

8)
Now lets see the files and folders created by SNAPPy during the subtyping process::

    ls

9)
As you can see the folders 'aligned'; 'blast'; 'trees' and the files 'subtype_results.csv' and 'report_subtype_results.csv' were created. The file 'subtype_results.csv' has the results that we wanted from running these task, the file 'report_subtype_results.csv' has a more extensive description of the tool outputs (more about this in the Tutorial :ref:`tut3`). The folders 'aligned', 'blast', and 'trees' contain the intermediate files created by SNAPPy during the subtyping process ( more about this in the section :ref:`how_it_works`). Now lets look at the file 'subtype_results.csv' (Note: You can use any other text editor or spreadsheet visualization tool.)::

    cat subtype_results.csv



10)
As you can see in this csv file there is the header of each fasta in the input followed by the result from SNAPPy subtyping. The 'id' numbers refers to the internal identifier used during the pipeline and links to the intermediate files in the 'aligned', 'blast', and 'trees' folders.


11)
Lets try to run the exactly same task again::

    snakemake subtype_all


Expected result::

    Building DAG of jobs...
    Nothing to be done.

12)
Nothing was done because the output that we requested was already built! This is one of the great advantages of using a pipeline software like `Snakemake <https://snakemake.readthedocs.io/en/stable/index.html>`_ , it goes top down looking for the requested file and the files needed to create it. If it is already there nothing needs to be done.

13)
Now lets use the SNAPPy rule thst allows us to clean all the outputs from a previous run. Attention!! If you have results that you want to keep change their name or move them to another folder before running the clean-up command::

    snakemake delete_all_outputs
    ls

14)
We are now back where we started without any output built. Lets run the pipeline but this time lets use more computational resources, namely four cpu threads::

    snakemake subtype_all --cores 4

15)
As you probably noticed this time the same process took a lot less time to run, that's because SNAPPy leverages the `Snakemake <https://snakemake.readthedocs.io/en/stable/index.html>`_ capabilities of parallelizing tasks. This allows SNAPPy to be extremely scalable. For instance if you have accesses to a n core cpu in theory you can use all of them to do subtyping with SNAPPy in one single task.

16)
That's it for this tutorial! If you now want to use SNAPPy on your own date don't forget to clean the outputs created during this tutorial and adjust the content of the input folder.

 .. _tut2:

2: Alignments with SNAPPy
^^^^^^^^^^^^^^^^^^^^^^^^^

Disclaimer: This tutorial and Tutorial :ref:`tut1` share the first four steps.

Al thought SNAPPy was built for HIV-1 subtyping one of its intermediary tasks is alignment to the reference genome (`HXB2 (K03455) <https://www.hiv.lanl.gov/components/sequence/HIV/asearch/query_one.comp?se_id=K03455>`_). Since SNAPPy is built based on `Snakemake <https://snakemake.readthedocs.io/en/stable/index.html>`_ we can call intermediary tasks, such as alignment, without running the entire pipeline. Making SNAPPy extremely useful for performing HIV-1 alignments.
To use SNAPPy the terminal always needs to be executed from the snappy folder, and the conda environment is active::

    conda actiavte snappy

We highly encourage to run the tests before using the pipeline just to make sure the tasks will run as intended and the pipeline was not modified or damaged::

    py.test

1)
Look at the content of the snappy folder::

    ls

Expected result::

    config.yaml  data  environment.yaml  input  LICENSE.txt  README.md  scripts  Snakefile  snappy_installer_x86_64.sh  test  test_pipeline.py

 
2)
In this tutorial we will focus on the folder 'input' and the  file 'config.yaml'. For exploration of the remaining content please look at the :ref:`how_it_works` section. Look inside the folder input::

    ls input

Expected result::

    test_msa.fasta

3)
As you can see SNAPPy has a file inside the folder 'input'. This is a test file used in this tutorial. For a successful run only one file can be inside the input folder. Therefore, if you want to run an analysis with your data you need to first delete the default file inside the input folder.

4)
Now lets look at the 'config.yaml' file:´

    cat config.yaml


Expected result::

    genomic_region:
      'GAG-POL-ENV'

5)
As you can see the genomic region selected for the SNAPPy's analysis is 'GAG-POL-ENV' which basically consists the HIV-1 regions named GAG, POL, and ENV concatenated ( according to the `HXB2 (K03455) <https://www.hiv.lanl.gov/components/sequence/HIV/asearch/query_one.comp?se_id=K03455>`_ reference genome). Further on this tutorial we will use different genomic regions.

6)
In order to perform the alignment of the sequences in the folder 'input' for the region specified in 'config.yaml' simple type::

    snakemake align_all

Expected last line of the result::

    5 of 5 steps (100%) done


7)
As you can see SNAPPy started executing the tasks needed (at green) to obtain the requested multiple sequence alignment (MSA). Lets see which files and folders SNAPPy created::

    ls

8)
The folder 'aligned' and the file 'all_aligned.fasta' were created. The folder 'aligned' contains intermediate files used to create the final MSA. You can now use a text editor or your favorite FASTA file reader (for instance `AliView <https://ormbunkar.se/aliview/>`_) to look at the 'all_aligned.fasta' file.As you can see it contains a lot of gaps ('-') because the aligned sequences only contained information for the GAG region and we requested an alignment to the HXB2 reference genome for the GAG,POL, and ENV regions ( as specified in the 'config.yaml file'). The produced sequences are of length 6918 nucleotides.

9)
Now lets save the obtained alignment to a new folder called ‘safe_outputs’ and use a SNAPPy rule to clean all the outputs previously created::

    mkdir safe_outputs
    cp all_aligned.fasta  safe_outputs/msa_gag_pol_env.fasta
    snakemake delete_all_outputs

10) 
The SNAPPy rule 'delete_all_outputs' is extremely useful to quickly delete files from previous runs but make sure that if you want to save outputs that you want to keep before running this rule (as we did above).

11)
Lets modify the 'config.yaml' file to obtain an alignment only for the GAG region. Open The 'config.yaml' in your favorite text editor and edit it so it looks like this::

    genomic_region:
      'GAG'

12)
Now lets ask SNAPPy to align the sequences in the input folder::

    snakemake align_all

13)
The same outputs as before were created. If we now evaluate the 'all_aligned.fasta' file we can see it has far less gaps ('-'). As stated before the inputs sequences only contained information for the GAG region, and these outputs (of length 1503 nucleotides) only have that said region.

14)
Fell free to test SNAPPy to create MSA for other HIV-1 sequences or using other genomic regions. Don't forget to always clean and save the outputs from previous runs if you want to keep them. If you are planning on exploring with different genomic regions don't forget to edit the 'config.yaml file'. The implemented genomic regions in SNAPPy are::

    'GAG', 'PR', 'RT', 'PR-RT',
    'INT', 'POL', 'ENV',
    'GAG-POL-ENV' 
    

15)
That's it for this tutorial! Don't forget that if you plan on using SNAPPy for subtyping the 'config.yaml' file always needs to indicate the option 'GAG-POL-ENV'.

 .. _tut3:

3: Results Analysis
^^^^^^^^^^^^^^^^^^^

In this tutorial we will give a more in-depth look at the outputs created by SNAPPy in the subtyping process. This tutorial starts after Tuturial :ref:`tut1`, and uses the outputs created in that tutorial. If you have not run Tutorial 1 yet or no longer have its outputs in the folder please do so before the next steps.



1)
Look at the content of the snappy folder::

    ls

Expected result::

    aligned  blast  config.yaml  data  environment.yaml  input  LICENSE.txt  README.md  scripts  Snakefile  snappy_installer_x86_64.sh  report_subtype_results.csv subtype_results.csv  test  test_pipeline.py  trees 

 
2)
In this tutorial we will focus on the files 'report_subtype_results.csv' and 'subtype_results.csv'. To read and edit them fell free to use your favorite text editor or spreadsheet reader. Lets open the 'subtype_results.csv' file.

3)
This is an extremely simple file with only tree columns: 'id', 'name', 'result'.

3.1)
The 'id' field in only important if you want to evaluate by yourself the intermediate files created by SNAPPy in the 'aligned', 'blast' and 'trees' folders. For instance the files refering to the FASTA with the header 'test01' will be named '0' like: 'aligned/0.fasta', 'aligned/aligned_0.fasta', 'blast/blast_0.txt', 'blast/recblast_0.txt', 'trees/all_0.nwk', 'trees/pure_0.nwk', and 'trees/recomb_0.nwk'. If you want to know more on why and how those files were created please see the :ref:`how_it_works` section.  

3.2)
The 'name' field corresponds to the headers found in the HIV-1 sequences in the file inside the 'input' folder. This will be the field that allows the user to cross the SNAPPy outputs with the user nomenclature.

3.3)
The 'result' field only contains the output produced by SNAPPy regathering that FASTA sequence subtype. No information is displayed regathering the analysis by BLAST or phylogenetic inference or how the decision was made. That information is in the file 'report_subtype_results.csv'.

4)
Now lets open the report_subtype_results.csv. This output has 12 columns named: 'id', 'name', 'result', 'recomb_result', 'node_all_refs', 's_node_all_refs', 'node_pure_refs', 's_node_pure_refs', 'node_recomb_refs', 's_node_recomb_refs', 'closser_ref', and 'rule'. The first three are exactly the same as for the subtype_results.csv file ( explained it points 3.1 to 3.3). The remaining will be described in the following topics.

4.1)
The field 'recomb_result' referes to an output obtained in a sliding with multiple BLASTs of the input. What this means is that the input was sliced multiple times and each slice served as a BLAST input to a database containing HIV-1 reference sequences. This test was mainly done looking for evidence of recombination in the target sequence. If you want to read more about the sliding window applied or the reference sequences used please read the :ref:`how_it_works` section. 

4.2)
The fields 'node_all_refs', 'node_pure_refs', and 'node_recomb_refs' correspond to the output of the phylogenetic inference using FastTree. These fields demonstrate if the target sequence was in a monophyletic clade with a group of HIV-1 reference sequences of only one subtype/circulation recombinant form (CRF). As the name indicate the 'node_all_refs' was inferred from a phylogenetic tree with potentially references from subtypes and CRFS, for the 'node_pure_refs' only subtype references were present, and for the 'node_recomb_refs' only CRFs references were present. If you want to know more about the parameters used in these phylogenetic inferences and/or the references used please go to the section :ref:`how_it_works`.

4.3)
The fields 's_node_all_refs', 's_node_pure_refs', and 's_node_recomb_refs' contain the support values for the monophyletic nodes where the criteria explain in point 4.2 are meet. To obtain these support values the Shimodaira-Hasegawa test as implemented in FastTree was used.

4.4)
The field 'closser_ref' shows the subtype or CRF of the reference sequence that showed to be closer to the target sequence in a BLAST analysis with all HIV-1 reference sequences used. If you want to know more about this BLAST or the reference sequences used please go to the section :ref:`how_it_works`.

4.5)
The field 'rule' is merely informative and show which SNAPPy 'rule' was used to make the decision about the subtyping output in the result field based on the other fields. If you want to know more about these rules please go to the section :ref:`how_it_works`.

5)
HIV-1 subtyping is a somewhat an ambiguous process. We believe that simplifying the 'subtype_results.csv' file allows users to quickly use SNAPPy, while providing the 'report_subtype_results.csv' file allows the user to observe the intermediate results created by SNAPPy and the decisions made. Please keep in mind that there will be cases harder to subtype that others, but you can always come back to the report_subtype_results.csv file and understand why SNAPPy outputted a given result.

6)
After this tutorial we believe that you are equipped with the knowledge to use SNAPPy and completely understand its outputs.

