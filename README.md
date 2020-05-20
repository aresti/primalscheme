# primalscheme

![Build](https://github.com/aresti/primalscheme/workflows/Build/badge.svg)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
---
A primer3 wrapper for designing multiplex primer schemes.

### requirements
* Linux, macOS or Windows
* Python 3.5+


### installation from PyPI using pip
Create & activate a virutal env (recommended):
```
python3 -m venv /path/to/primal-venv
source /path/to/primal-venv/bin/activate
```
or on Windows:
```
python -m venv c:\path\to\primal-venv
c:\path\to\primal-venv\Scripts\activate.bat
```

Install primalscheme with pip:
```
pip install primalscheme
```

### installation from source
```
git clone https://github.com/aresti/primalscheme.git primalscheme
cd primalscheme

python3 -m venv venv
source venv/bin/activate

pip install .
```

Note: if you'd like your installation to be editable, use:
```
pip install flit
flit install --symlink
```

### usage
Download our CHIKV_demo.fa file from <a href="https://raw.githubusercontent.com/aresti/primalscheme/master/tests/inputs/CHIKV_demo.fa">tests/inputs/CHIKV_demo.fa</a>

```
primalscheme multiplex path/to/CHIKV_demo.fa
```
or
```
python -m primalscheme multiplex path/to/CHIKV_demo.fa
```

### about
Primal Scheme is a tool for designing multiplex PCR primers for generating tiling amplicons. It was developed for sequencing large numbers of viral isolates of known lineages e.g. outbreak strains. It requires only two PCR reactions to generate the products which should cover at least 90% of the target region without optimisation. Full coverage is possible by fine tuning of primer concentration or exchanging poor performers.

The primal scheme software is a wrapper for primer3 which is used to generate candidate primers from the primary reference (the first reference in the FASTA file). The minimum requirements are a single reference in FASTA format but if additional references are supplied it then aligns the primers to all reference genomes before selecting the pair with the highest mean identity. The input FASTA file should contain complete genomes representative of the lineages you would expect to find in your samples. The first genome in the fasta file is used to generate the candidates so it is the most important, see guidelines below.

Most parameters for running primer3 e.g. Tm are not directly exposed in the CLI as we believe they are necessary for successful multiplex PCR. We strongly recommend using the PCR conditions outlined in the <a href="http://www.nature.com/nprot/journal/v12/n6/full/nprot.2017.066.html">Nature Protocols paper</a>. However, most values can be modified in `config.json` if you install the package from source with `flit install -s` (see above). We have sequenced viral genomes of up to 20 kb in size with amplicon sizes from 300 to 1000 bp. You can adjust the desired amplicon size (min and max) to suit your sample type or sequencing platform. For example, using an amplicon size of 400 (min 380, max 420) suits both the MiSeq 2x250 run configuration or MinION. You can also adjust the desired overlap parameter to set the size of sequence shared by overlapping amplicons, although it is set to 0 by default. If the desired overlap cannot be found it will fall back to the largest available.

The primers required will be output into a .TSV file containing additional information such as size, %GC and Tm. Reference genomes and primer positions are written to .FASTA and .BED file which are required for primer trimming downstream (see: https://github.com/artic-network/fieldbioinformatics). Other files useful for assessing output are the .LOG (use --debug for verbose logging), .SVG/.PDF show a diagramatic representation of the primers scheme and .PICKLE which contains the primalscheme objects including alternative primer options. 

Guidelines for designing a scheme:

1. Download complete genomes from GenBank and put into a single FASTA file
2. Align all sequences using Clustal Omega
3. Check start and ends of the alignment and remove any sequences that are not full size
4. Click the 'Phylogenetic Tree' tab to see how many major clades there are
5. Click the 'Result Summary' tab then 'Percent Identity Matrix' link to see the identity matrix of all sequences
6. Check that the maximum sequence divergence is not more than 10%, it may be necessary to group the more similar genomes and separate them into multiple schemes
7. Remove duplicate or sequences with >99% identity to other genomes in the file as these will bias the primer selection
8. Ensure the first genome in the file is the most recent, complete genome for a given geographical region/lineage
9. Run primal scheme on the FASTA file