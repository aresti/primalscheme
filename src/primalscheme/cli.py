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

import argparse
import logging
import sys

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from progress.bar import ShadyBar

from primalscheme import __version__ as version, config
from primalscheme.region import NoSuitablePrimersError
from primalscheme.reporting import MultiplexReporter, ProgressTracker

logger = logging.getLogger("primalscheme")


def main():
    """
    Entry point for primalscheme.
    Process args, set up logging and call the command function.
    """

    # Parse args, update config
    try:
        args = parse_arguments(sys.argv[1:])
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    # Validate output path
    try:
        outpath = get_output_path(args.output_path, force=args.force)
    except IOError as e:
        print(f"Error: {e}")
        sys.exit(1)

    # Setup logging
    setup_logging(outpath, debug=args.debug, prefix=args.prefix)

    for arg in vars(args):
        logger.debug("{}: {}".format(arg, str(vars(args)[arg])))

    # Run subcommand
    args.func(args, outpath)


def multiplex(args, outpath):
    """Multipex scheme command."""

    # Process FASTA input
    min_ref_size = args.amplicon_size_max * 2
    sort = not args.no_sort
    try:
        references = process_fasta(args.fasta, min_ref_size=min_ref_size, sort=sort)
    except ValueError as e:
        logger.error(f"Error: {e}")
        sys.exit(2)

    # Progress bar
    progress_bar = ProgressBar()

    # Create scheme
    try:
        scheme = MultiplexReporter(
            outpath,
            references,
            prefix=args.prefix,
            amplicon_size_min=args.amplicon_size_min,
            amplicon_size_max=args.amplicon_size_max,
            target_overlap=args.target_overlap,
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
    if args.debug:
        scheme.write_pickle()
    sys.exit(0)


def process_fasta(file_path, min_ref_size=None, sort=True):
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

    # Sort references by length to maximum scheme coverage
    if sort:
        references.sort(key=len, reverse=True)
    return references


def setup_logging(output_path, debug=False, prefix="primalscheme"):
    """Setup logging output and verbosity."""

    logger.setLevel(logging.DEBUG if debug else logging.INFO)

    # File handler
    log_filepath = output_path / f"{prefix}.log"
    fh = logging.FileHandler(log_filepath)
    fh.setLevel(logging.DEBUG)
    fh_formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    fh.setFormatter(fh_formatter)
    logger.addHandler(fh)

    # Stream handler STDOUT
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.INFO)
    sh_formatter = logging.Formatter("%(message)s")
    sh.setFormatter(sh_formatter)
    logger.addHandler(sh)

    logger.info(f"Writing log to {log_filepath}")


def get_output_path(output_path, force=False):
    """
    Check for an existing output dir, require --force to overwrite.
    Create dir, return path object.
    """
    path = Path(output_path)

    if path.exists() and not force:
        raise IOError("Directory exists add --force to overwrite")
    elif path.exists() and not path.is_dir():
        raise IOError("The output path is not a directory.")

    path.mkdir(exist_ok=True)
    return path


def parse_arguments(args):
    """Parse command line arguments."""
    # Setup parsers
    parser = argparse.ArgumentParser(
        prog="primalscheme",
        description=("A primer3 wrapper for designing multiplex primer schemes."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparsers = parser.add_subparsers(title="[subcommands]", dest="command")
    subparsers.required = True

    parser_scheme = subparsers.add_parser("multiplex", help="Multiplex PCR scheme")

    # Set func
    parser_scheme.set_defaults(func=multiplex)

    # Add arguments to parser
    parser_scheme.add_argument("fasta", help="FASTA file")
    parser_scheme.add_argument(
        "--prefix",
        default=config.PREFIX,
        help="Prefix used for primer names and output files " "(default: %(default)s)",
    )
    parser_scheme.add_argument(
        "--amplicon-size-min",
        type=positive_int,
        default=config.AMPLICON_SIZE_MIN,
        help="Minimum amplicon size (default: %(default)i)",
    )
    parser_scheme.add_argument(
        "--amplicon-size-max",
        type=positive_int,
        default=config.AMPLICON_SIZE_MAX,
        help="Maximum amplicon size (default: %(default)i)",
    )
    parser_scheme.add_argument(
        "--output-path",
        default=config.OUTPUT_PATH,
        help="Output directory (default: %(default)s)",
    )
    parser_scheme.add_argument(
        "--target-overlap",
        default=config.TARGET_OVERLAP,
        type=positive_int,
        help="Target overlap size (default: %(default)i)",
    )
    parser_scheme.add_argument(
        "--step-distance",
        default=config.STEP_DISTANCE,
        type=positive_int,
        help="Distance to step between find attempts (default: %(default)i)",
    )
    parser_scheme.add_argument("--debug", action="store_true", help="Verbose logging")
    parser_scheme.add_argument(
        "--no-sort",
        action="store_true",
        help="Don't sort input FASTA by length (will use first reference in BED file",
    )
    parser_scheme.add_argument(
        "--force",
        action="store_true",
        help="Force output to an existing directory and overwrite output " "files",
    )
    parser.add_argument(
        "-V", "--version", action="version", version=f"%(prog)s {version}"
    )

    # Generate args
    parsed = parser.parse_args(args)

    # Post parsing checks
    if parsed.amplicon_size_max < parsed.amplicon_size_min:
        raise ValueError("--amplicon-size-min cannot exceed --amplicon-size-max")

    return parsed


def positive_int(string):
    """Validate command line args as being of type: positive int."""
    value = int(string)
    if value < 0:
        raise argparse.ArgumentTypeError("positive integer required.")
    return value


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
    main()
