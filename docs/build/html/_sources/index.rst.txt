.. snappy-docs documentation master file, created by
   sphinx-quickstart on Thu Jun 13 10:52:37 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   

Welcome to snappy's documentation!
==================================



 
   
SNAPPy's code can be found `here <https://github.com/PMMAraujo/snappy>`_ 



**What is SNAPPy**


SNAPPy is a Snakemake pipeline for HIV-1 subtyping by phylogenetic pairing. It was build to allow locaal high-throughput HIV-1 subtyping while being resource efficient and scalable. This pipeline was constructed using `Snakemake <https://snakemake.readthedocs.io/en/stable/index.html>`_ , and it uses `MAFFT <https://mafft.cbrc.jp/alignment/software/>`_ for multiple sequence alignments, `BLAST <https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs>`_ for similarirty queries, `IQ-TREE <http://www.iqtree.org/>`_ for phylogenetic inference, and several `Biopython <https://biopython.org/>`_ modules and the `ETE toolkit <http://etetoolkit.org/>`_ for data parsing an analysis. For in-depth information on how the tool works please visit the :ref:`how_it_works` section. SNAPPy can be used in Linux based distributions (tested on Ubuntu 16 and 18), macOS 10 (tested of Mojave 10.14), and Windows 10 using Windows Subsystem for Linux (tested with Ubuntu 18 bash).


**Getting Started**


* Installation:

   * :ref:`quick_l`

   * :ref:`quick_m`

   * :ref:`point-by-point_l`

   * :ref:`point-by-point_m`

* Tutorials:

   * :ref:`tut1`

   * :ref:`tut2`

   * :ref:`tut3`

* Commands list:

   * :ref:`usage`


**Support**

 * The software is provided "as is", without warranty of any kind, we do our best to keep it updated and bug free
 *  :ref:`faq`
 * For additional information please contact: INFO (MAIL) 


**How to Cite**

text



.. toctree::
   :hidden:
   :maxdepth: 2
   :titlesonly:
   
   installation
   tutorials
   usage
   how_it_works
   faq
   cite
   license
 
