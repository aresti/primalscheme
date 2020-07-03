"""
PrimalScheme: a primer3 wrapper for designing multiplex primer schemes

Copyright (C) 2020 Joshua Quick and Andrew Smith
www.github.com/aresti/primalscheme

This module contains the CLI for PrimalScheme.
It is executed when the user runs 'primalscheme' after installation,
or 'python -m primalscheme'.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>
"""

import click
import logging
import sys

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from progress.bar import ShadyBar

from primalscheme import __version__ as version, config
from primalscheme.primer import calc_gc
from primalscheme.region import NoSuitablePrimersError
from primalscheme.reporting import MultiplexReporter, ProgressTracker

logger = logging.getLogger("primalscheme")


CLI_CONTEXT = dict(auto_envvar_prefix="PRIMAL", help_option_names=["-h", "--help"],)


@click.group(context_settings=CLI_CONTEXT)
@click.version_option(version, "--version", "-V")
def cli():
    """a tool for designing primer panels for multiplex PCR."""
    pass


@cli.command()
@click.argument("fasta", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--amplicon-size",
    "-a",
    multiple=True,
    type=click.IntRange(90),
    help=(
        "Amplicon size target. Pass twice to set an exact range, otherwise expect "
        f"+/- {config.SIZE_RANGE_AUTO/2 * 100}%."
    ),
    metavar="<int>",
    default=[config.AMPLICON_SIZE_MIN, config.AMPLICON_SIZE_MAX],
    show_default=True,
)
@click.option(
    "--outpath",
    "-o",
    type=click.Path(file_okay=False, writable=True),
    help="Path to output directory.",
    metavar="<dir>",
    default=config.OUTPUT_PATH,
    show_default=True,
)
@click.option(
    "--name",
    "-n",
    type=click.STRING,
    help="Prefix name for your outputs.",
    metavar="<str>",
    default=config.PREFIX,
    show_default=True,
)
@click.option(
    "--target-overlap",
    "-t",
    type=click.IntRange(0),
    help="Target insert overlap size.",
    metavar="<int>",
    default=config.TARGET_OVERLAP,
    show_default=True,
)
@click.option("--debug/--no-debug", "-d", help="Set log level DEBUG.", default=False)
@click.option(
    "--force/--no-force",
    "-f",
    help="Force output to an existing directory, overwrite files.",
    default=False,
)
@click.option(
    "--pinned/--no-pinned",
    "-p",
    help="Only consider primers from the first reference.",
    default=False,
)
@click.option(
    "--high-gc/--no-high-hc",
    "-g",
    help="Use config suitable for high-GC sequences.",
    default=False,
)
def multiplex(
    fasta, amplicon_size, outpath, name, target_overlap, debug, force, pinned, high_gc
):
    """Design a multiplex PCR scheme."""
    # Handle args
    if len(amplicon_size) == 1:
        target_amplicon_size = amplicon_size[0]
        half_range = int(target_amplicon_size * config.SIZE_RANGE_AUTO / 2)
        amplicon_size_min = target_amplicon_size - half_range
        amplicon_size_max = target_amplicon_size + half_range
    else:
        amplicon_size_min = min(amplicon_size[:2])
        amplicon_size_max = max(amplicon_size[:2])

    # Validate output path
    try:
        outpath = get_output_path(outpath, force=force)
    except IOError as e:
        click.echo(click.style(f"Error: {e}", fg="red",))
        sys.exit(1)

    # Setup logging
    setup_logging(outpath, debug=debug, prefix=name)

    # Process FASTA input
    try:
        references = process_fasta(fasta, min_ref_size=amplicon_size_max)
    except ValueError as e:
        logger.error(f"Error: {e}")
        sys.exit(2)

    # High GC warning
    if not high_gc and any(
        calc_gc(str(ref.seq)) > config.HIGH_GC_WARN_THRESHOLD for ref in references
    ):
        click.echo(
            click.style(
                "WARNING: High-GC reference detected. Consider using --high-gc mode.",
                fg="red",
            )
        )

    # Progress bar
    progress_bar = ProgressBar()

    # Log references
    ref_ids = [f" - {ref.id}" for ref in references]
    logger.info("\n".join(["References:"] + ref_ids))
    logger.info(f"Primary reference for coordinate system: {references[0].id}")
    logger.info(
        "Considering primers from "
        f"{'primary reference only' if pinned else 'all references'}"
    )

    # Create scheme
    try:
        scheme = MultiplexReporter(
            outpath,
            references,
            prefix=name,
            amplicon_size_min=amplicon_size_min,
            amplicon_size_max=amplicon_size_max,
            target_overlap=target_overlap,
            primary_only=pinned,
            high_gc=high_gc,
            progress_tracker=progress_bar,
        )
        scheme.design_scheme()
    except NoSuitablePrimersError:
        # Unable to find primers for at least one region
        logger.error("Error: Unable to find suitable primers")
        sys.exit(10)
    except Exception as e:
        # Unexpected error
        logger.error(f"Error: {e}")
        raise e

    # Write outputs
    logger.info(
        f"All done! Scheme created with {len(scheme.regions)} regions, "
        f"{scheme.gap_count } gap{'' if scheme.gap_count == 1 else 's'}, "
        f"{scheme.percent_coverage}% coverage"
    )
    scheme.write_default_outputs()
    if debug:
        scheme.write_pickle()
    sys.exit(0)


def process_fasta(file_path, min_ref_size=None):
    """Parse and validate the fasta file."""

    references = []
    records = SeqIO.parse(file_path, "fasta")  # may raise

    # Remove gaps
    for record in records:
        ref = SeqRecord(
            Seq(str(record.seq).replace("-", "").upper()),
            id=record.id,
            description=record.id,
        )
        references.append(ref)

    # Check for no references
    if not references:
        raise ValueError("The input FASTA file does not contain any valid references.")

    # Check for too short references
    if min_ref_size and any(len(ref) < min_ref_size for ref in references):
        raise ValueError(
            "One or more of your references is too short. Based on your target "
            f"amplicon size, the minimum reference size is {min_ref_size} nt."
        )

    # Check for too many references
    if len(references) > 100:
        raise ValueError("A maximum of 100 reference genomes is currently supported.")

    # Check for max difference in size between references
    primary_ref = references[0]
    primary_ref_len = len(primary_ref)

    if any(abs(len(r) - primary_ref_len) > 500 for r in references):
        raise ValueError(
            "One or more of your references is too different in size to "
            "the primary (first) reference. The maximum difference is "
            "500 nt."
        )

    # Check for a valid alphabet
    VALID_ALPHABET = "AGCTacgt"
    for r in references:
        if any(base not in VALID_ALPHABET for base in set(r.seq)):
            raise ValueError(
                "One or more of your fasta sequences contain invalid "
                f"nucleotide codes. The supported alphabet is '{VALID_ALPHABET}'. "
                "Ambiguity codes and gaps are not currently supported."
            )

    return references


def setup_logging(output_path, debug=False, prefix="primalscheme"):
    """Setup logging output and verbosity."""

    logger.setLevel(logging.DEBUG if debug else logging.INFO)

    # File handler
    log_filepath = output_path / f"{prefix}.log"
    fh = logging.FileHandler(log_filepath)
    fh.setLevel(logging.DEBUG)
    fh_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    fh.setFormatter(fh_formatter)
    logger.addHandler(fh)

    # Stream handler STDOUT
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.INFO)
    sh_formatter = logging.Formatter("%(message)s")
    sh.setFormatter(sh_formatter)
    logger.addHandler(sh)

    logger.info(f"Writing log to {log_filepath}")
    logger.debug(f"primalscheme version {version}")


def get_output_path(output_path, force=False):
    """
    Check for an existing output dir, require --force to overwrite.
    Create dir, return path object.
    """
    path = Path(output_path)

    if path.exists() and not force:
        raise IOError("Directory exists add --force to overwrite")

    path.mkdir(exist_ok=True)
    return path


class ProgressBar(ShadyBar, ProgressTracker):
    """Progress bar for terminal stdout."""

    suffix = "%(percent)d%% [%(index)d / %(max)d]"

    def __init__(self, *args, **kwargs):
        """Init ProgressBar."""
        self.__considered = 0
        self.__region_num = 1
        super().__init__(*args, **kwargs)

    def goto(self, val):
        """Update progress to val."""
        super().goto(val)

    @property
    def end(self):
        """End (max) progress value."""
        return self.max

    @end.setter
    def end(self, val):
        """Set end progress value."""
        self.max = val

    @property
    def considered(self):
        """Count of considered primers."""
        return self.__considered

    @considered.setter
    def considered(self, considered):
        """Set count of considered primers."""
        self.__considered = considered
        self.update_message()

    @property
    def region_num(self):
        """Current region num."""
        return self.__region_num

    @region_num.setter
    def region_num(self, region_num):
        """Set current region num."""
        self.__region_num = region_num
        self.update_message()

    def update_message(self):
        """Update progress bar prefix message."""
        self.message = f"Considered {self.considered} primers, region {self.region_num}"

    def interrupt(self):
        """Prepare to be interrupted by a log message."""
        if self.index:
            self.finish()


if __name__ == "__main__":
    cli()
