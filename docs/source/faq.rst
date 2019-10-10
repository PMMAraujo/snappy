.. _faq:

FAQ
===


**Q: "Why should I use SNAPPy?**

A: Do you want to do HIV-1 subtyping/alignment in large numbers? Do you have to process your HIV-1 sequences locally because ethic concerns or data privacy policy? Do you want  a tool for HIV-1 subtyping/alignment that scales with the amount or computational resources you have? Do you want a tool that gives simple outputs but also allows you to dive deep in intermediate files and how the decisions were made? If you answered was yes to any or all or the previous questions SNAPPy is a tool for you.

**Q: "Can SNAPPy run on a windows machine?**

A: We tested SNAPPy in windows 10 using the `Windows Subsystem for Linux <https://docs.microsoft.com/en-us/windows/wsl/about>`_ . It should permor as well as in a native UNIX based machine.


**Q: "Is there any graphical user interface to use SNAPPy?**

A: Currently, SNAPPy is a command base tool. However, we are trying to implement more user friendly ways to use SNAPPy without losing functionality.

**Q: "Why using BLAST in SNAPPy?**

A: We decided to use BLAST for the first version of SNAPPy because it is reliable, extensively tested and easy to implement. We are aware of more "modern" BLAST-like tools, but one needs to recognize that BLAST has decades of being "battle tested". This decision does not restrict the usage of other tools/software in future versions of SNAPPy. 

**Q: "I used SNAPPy to subtype 100k sequences and the outputs are occupying a lot of disk space (several dozens gigabytes), is it normal?**

A: Please keep in mind that SNAPPy is a complex pipeline and performs several analysis. Each of those generate outputs. Everything is multiplied by the number of inputs, even if each input occupies an almost negligible amount of disk space the sum of all of them will result in a large value. At the end of each run SNAPPy automatically deletes an hidden folder created by Snakemake named '.snakemake', this folder is extremely large and its deletion will save some space ( do not delete it if the pipeline is still being executed). 


**Q: "Why can I not put several FASTA files inside the input folder?**

A: This was a development decision. This may change in the future but for the current SNAPPy version we believe that is the best option for both the user and the pipeline. 
