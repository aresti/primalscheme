# primal
A primer3 wrapper for designing multiplex primer schemes

Primal is a python package for generating multiplex tiling amplicon schemes for viral sequencing. As opposed to focusing on algorithms to generate the candidate primers i.e. K-mer matching or de Bruijn graphs, Primal uses Primer3 to generate candidate primers from a single reference genome. It then scores the primers based on pairwise alignment to any number of additional reference genome before selecting the best candidates.  It was the method we used to generate the Zika virus scheme used for the ZiBRA project after little success using existing methods.


