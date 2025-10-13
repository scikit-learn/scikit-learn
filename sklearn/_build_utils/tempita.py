#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""
Improved Tempita template processor.
"""

import argparse
import os
import sys

from Cython import Tempita as tempita


def process_tempita(fromfile: str, outfile: str | None = None) -> None:
    """Render a Tempita template and write the result.

    Parameters
    ----------
    fromfile : str
        Path to the input .tp file.
    outfile : str | None
        Path to the output file. If None, writes to stdout.
    """
    try:
        with open(fromfile, "r", encoding="utf-8") as f:
            template_content = f.read()
    except FileNotFoundError:
        sys.exit(f"Input file not found: {fromfile}")

    # Render the template
    try:
        template = tempita.Template(template_content)
        content = template.substitute()
    except Exception as e:
        sys.exit(f"Template rendering failed: {e}")

    # Write to file or stdout
    if outfile:
        os.makedirs(os.path.dirname(outfile), exist_ok=True)
        with open(outfile, "w", encoding="utf-8") as f:
            f.write(content)
        print(f"Generated: {outfile}")
    else:
        print(content)


def main() -> None:
    parser = argparse.ArgumentParser(description="Process Tempita templates (.tp)")
    parser.add_argument("infile", type=str, help="Path to the input file (.tp)")
    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        help="Path to the output directory (default: print to stdout)",
    )
    parser.add_argument(
        "-i",
        "--ignore",
        type=str,
        help="Optional ignored input, useful for build dependencies",
    )
    args = parser.parse_args()

    infile = args.infile
    if not infile.endswith(".tp"):
        sys.exit(f"❌ Unexpected file extension: {infile}")

    if args.outdir:
        outdir_abs = os.path.abspath(args.outdir)
        os.makedirs(outdir_abs, exist_ok=True)
        outfile_name = os.path.splitext(os.path.basename(infile))[0]
        outfile = os.path.join(outdir_abs, outfile_name)
    else:
        outfile = None

    process_tempita(infile, outfile)


if __name__ == "__main__":
    main()
