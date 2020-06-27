# primalscheme

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Build](https://github.com/aresti/primalscheme/workflows/Build/badge.svg)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=aresti_primalscheme&metric=alert_status)](https://sonarcloud.io/dashboard?id=aresti_primalscheme)

---

## about

primalscheme is a tool for designing primer panels for multiplex PCR. It uses a greedy algorithm to find primers for tiling amplicon generation for multiple reference genomes. It works best on viral isolates of known lineages e.g. outbreak strains.

## web interface

a web interface is available at [primalscheme.com](https://primalscheme.com/), although it does not currently support all available options.

## requirements

* Linux, macOS or Windows
* Python 3.6+

## installation from PyPI using pip

Create & activate a virtual environment (recommended):

```shell
python3 -m venv /path/to/primal-venv
source /path/to/primal-venv/bin/activate
```

or on Windows:

```shell
python -m venv c:\path\to\primal-venv
c:\path\to\primal-venv\Scripts\activate.bat
```

Install primalscheme with pip:

```shell
pip install primalscheme
```

## installation from source

```shell
git clone https://github.com/aresti/primalscheme.git primalscheme
cd primalscheme
python3 -m venv venv
source venv/bin/activate
pip install .
```

Note: if you'd like your installation to be editable, use:

```shell
pip install flit
flit install --pth-file
```

(`--pth-file` vs `--symlink`: maintains module visibility for IDE debuggers, and works on Windows.)

## command line interface

To design a multiplex scheme run `primalscheme multiplex <FASTA>`.

If you're looking for test inputs you can try [CHIKV_demo.fa](tests/inputs/CHIKV_demo.fa), [Ebov-10-Pan.fasta](tests/inputs/Ebov-10-Pan.fasta) or [nCov-2019.fasta](tests/inputs/nCov-2019.fasta)

### optional arguments

* `--amplicon-size, -a`: Amplicon size target. Pass twice to set an exact range (e.g, -a 220 -a 250), otherwise expect +/- 5% (default 380, 420).
* `--outpath, -o`: Path to output directory (default: ./output).
* `--name, -n`: Prefix name for your outputs (default: scheme).
* `--target-overlap, -t`: Target insert overlap size (default: 0).
* `--debug, -d`: Set log level DEBUG.
* `--force, -f`: Force output to an existing directory, overwrite files.
* `--pinned, -p`: Only consider primers from the first reference.
* `--help, -h`: Show help message and exit.

### modifying config

Many parameters are not directly exposed in the CLI as we believe they are necessary for successful multiplex PCR. However, most values can be modified in `config.py` if you install the package from source in editable mode with `flit install --pth-file` (see above).

### output

* `{output_dir}/{prefix}.reference.fasta` - all input references
* `{output_dir}/{prefix}.primer.bed` - coordinates of primer positions
* `{output_dir}/{prefix}.insert.bed` - coordinates of trimmed amplicons
* `{output_dir}/{prefix}.primer.tsv` - primer sequences and information
* `{output_dir}/{prefix}.plot.pdf` - diagrammatic overview of scheme (PDF)
* `{output_dir}/{prefix}.plot.svg` - diagrammatic overview of scheme (SVG)
* `{output_dir}/{prefix}.report.json` - run report (reference ids, regions, gaps, coverage)
* `{output_dir}/{prefix}.pickle` - pickled Python objects (if --debug)
* `{output_dir}/{prefix}.log` - run logs (more detail if --debug)

## how to design a scheme

1. Download complete genomes from GenBank and concatenate into a single FASTA file.
    * e.g. `cat "sequence (1).fasta" "sequence (2).fasta" > input.fasta`
    * primalscheme does not currently support inputs with ambiguity codes.
2. Align all sequences using [Clustal Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/).
3. The genomes should be complete and roughly the same length.
    * Primers are only designed within the limits of the alignment.
    * The **first reference** will be used for the coordinate system in the output scheme.
4. Check the 'Result Summary' tab then the 'Percent Identity Matrix' link to check the identity matrix of all sequences.
    * Check that the maximum sequence divergence is not more than 5%.
    * If your sequences are more divergent, it may be necessary to separate them into multiple independent schemes; the 'Phylogenetic Tree' tab may be useful for this.
5. Remove sequences with 99-100% identity to other genomes in the file as these will bias the primer selection
6. Run primalscheme on your FASTA file.
    * Optionally, set amplicon min/max length as required by your sequencing technology.

## protocol

We strongly recommend using the PCR conditions outlined in the [Nature Protocols paper](http://www.nature.com/nprot/journal/v12/n6/full/nprot.2017.066.html).

## description of the algorithm

* Set the first reference as primary, to be used for the coordinate system.
* Open a region window of width AMPLICON_SIZE_MAX onto the primary reference.
* Align the flanks of this reference slice to all other references using parasail.
* Take all aligned flanks and digest them into K-mers (primers) of length PRIMER_SIZE_MIN to PRIMER_SIZE_MAX.
* Hard filter the primers to eliminate those with unacceptable GC, Tm, homodimer and homopolymer length.
* Filter out any primers with > PRIMER_MAX_MISMATCHES against any reference.
* Sort left and right primers according to base penalty plus weighted mismatch score.
* Check for heterodimers in same pool using primer3-py; select top scoring, non-interacting pair.
* Open new window position set to maintain overlap if region is successful.
* Repeat routine until end of genome is reached.
* If during the above routine there is a) failed alignment b) no primers after filtering then:
  * Move window position left by STEP_DISTANCE before retrying.
  * Continue stepping left until a hard-left limit is reached (start of genome, primer collision or no progress).
  * Reset window position to initial value.
  * Move window position right (to open a gap or reduce overlap) before retrying.
  * Exclude secondary reference on repeated flank alignment failure.

## citing

If you use primalscheme please cite:

```Quick J et al. Multiplex PCR method for MinION and Illumina sequencing of Zika and other virus genomes directly from clinical samples. Nat Protoc. 2017 Jun;12(6):1261-1276.```
