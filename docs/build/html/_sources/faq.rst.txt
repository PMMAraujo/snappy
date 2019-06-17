 .. _faq:

FAQ
===


**Q: "Why should I use SNAPPy?**

A: Do you want to do HIV-1 subtyping/alignment in large numbers? Do you have to process your HIV-1 sequences locally because ethic concerns or data privacy policy? Do you want  a tool for HIV-1 subtyping/alignment that scales with the amount or computational resources you have? Do you want a tool that gives simple outputs but also allows you to dive deep in intermediate files and how the decisions were made? If you answered was yes to any or all or the previous questions SNAPPy is a tool for you.

**Q: "Can SNAPPy run on a windows machine?**

A: In theory it is possible to run SNAPPy on a windows machine using the `Windows Subsystem for Linux <https://docs.microsoft.com/en-us/windows/wsl/about>`_ , however at the moment we do not provide any support for this and if you decide to do so do it at your own risk. The main reason SNAPPy was build for Linux based systems is that is has dependencies that only work on these operating systems. We are considering options to make SNAPPy operative systems agnostic in future releases.


**Q: "Is there any graphical user interface to use SNAPPy?**

A: Currently, SNAPPy is a command base tool. However, we are trying to implement more user friendly ways to use SNAPPy without losing functionality.

**Q: "Why using BLAST in SNAPPy?**

A: We decided to use BLAST for the first version of SNAPPy because it is reliable, extensively tested and easy to implement. We are aware of more "modern" BLAST-like tools, but one needs to recognize that BLAST has decades of being "battle tested". This decision does not restrict the usage of other tools/software in future versions of SNAPPy. 

**Q: "Why using FastTree in SNAPPy?**

A: After testing several phylogenetic tools we acknowledged that FasTree had the best ratio between phylogenetic inference quality and computational time. Moreover, FastTree has all the characteristics needed for the  phylogenetic inference software for this version of SNAPPy. This decision does not restrict the usage of other tools/software in future versions of SNAPPy. 

**Q: "I used SNAPPy to subtype 100k sequences and the outputs are occupying a lot of disk space (several dozens gigabytes), is it normal?**

A: Please keep in mind that SNAPPy is a complex pipeline and performs several analysis. Each of those generates outputs. Everything is multiplied by the number of inputs, even if each input occupies an almost negligible amount of disk space the sum of all of them will result in a large value. At the end of each run SNAPPy automatically deletes  an hidden folder created by Snakemake named '.snakemake', this folder is extremely large and its deletion will save some space ( do not delete if the pipeline is still being executed). 


**Q: "Why can I not put several FASTA files inside the input folder?**

A: This was a development decision. This may change in the future but for the current SNAPPy version we believe thats the best option for both the user and the pipeline. 
