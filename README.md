# primalscheme
A primer3 wrapper for designing multiplex primer schemes.

# requirements
* macOS or Linux
* Python 3.6
* C++ compiler (for the porechop dependency)
    - If you're using GCC, version 4.9.1 or later is required (check with g++ --version)

Not tested on Windows, but it should be possible to install.


# installation
```
git clone https://github.com/aresti/primalscheme.git
cd primalscheme
virtualenv -p python3 venv
source venv/bin/activate
pip install .
```

Or, from within an existing an environment simply:
```
pip3 install git+https://github.com/aresti/primalscheme.git
```

# usage
```
primalscheme multiplex tests/inputs/CHIKV_demo.fa
```

# about
Primal Scheme is a tool for designing multiplex PCR primers for generating tiling amplicons. It was developed for sequencing large numbers of viral isolates of known lineages e.g. outbreak strains. It requires only two PCR reactions to generate the products which should cover at least 80% of the target region without optimisation. Higher coverage is possible with balancing of individial primer concentrations.

The primal scheme software is a wrapper for primer3 which is used to generate candidate primers from the primary reference (the first reference in the FASTA file). It then scores the primers based on pairwise alignment to all additional reference genomes before selecting the most universal. The input FASTA file should contain complete genomes representative of the lineages you would expect to find in your samples. The first genome in the fasta file is used to generate the candidates so it is the most important, see guidelines below.

Parameters for running primer3 e.g. length and Tm are hardcoded as we believe they are necessary for successful multiplex PCR. We strongly recommend using the PCR conditions outlined in the <a href="http://www.nature.com/nprot/journal/v12/n6/full/nprot.2017.066.html">Nature Protocols paper</a>. We have sequenced viral genomes of up to 20 kb in length with amplicon lengths from 300 to 1000 bp. You can adjust the desired product length to suit your sample type or sequencing platform. For example, using an amplicon length of 400 suits both the MiSeq 2x250 run configuration or MinION. You can also adjust the overlap parameter to set the length of sequence shared by overlapping amplicons.

Guidelines for designing a scheme:

1. Download complete genomes from GenBank and put into a single FASTA file
2. Align all sequences using Clustal Omega
3. Check start and ends of the alignment and remove any sequences that are not full length
4. Click the 'Phylogenetic Tree' tab to see how many major clades there are
5. Click the 'Result Summary' tab then 'Percent Identity Matrix' link to see the identity matrix of all sequences
6. Check that the maximum sequence divergence is not more than 10%, it may be necessary to group the more similar genomes and make separate schemes for them
7. Remove duplicate or sequences with >99% identity to other genomes in the file as these will bias the primer selection
8. Ensure the first genome in the file is the most recent, complete genome for a given geographical region/lineage
9. Run primal scheme on the FASTA file
