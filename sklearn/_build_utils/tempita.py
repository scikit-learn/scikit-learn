#!/usr/bin/env python3

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import argparse
import os

from Cython import Tempita as tempita

# XXX: If this import ever fails (does it really?), vendor either
# cython.tempita or numpy/npy_tempita.


def process_tempita(fromfile, outfile=None):
    """Process tempita templated file and write out the result.

    The template file is expected to end in `.c.tp` or `.pyx.tp`:
    E.g. processing `template.c.in` generates `template.c`.

    """
    with open(fromfile, "r", encoding="utf-8") as f:
        template_content = f.read()

    template = tempita.Template(template_content)
    content = template.substitute()

    with open(outfile, "w", encoding="utf-8") as f:
        f.write(content)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", type=str, help="Path to the input file")
    parser.add_argument("-o", "--outdir", type=str, help="Path to the output directory")
    parser.add_argument(
        "-i",
        "--ignore",
        type=str,
        help=(
            "An ignored input - may be useful to add a "
            "dependency between custom targets"
        ),
    )
    args = parser.parse_args()

    if not args.infile.endswith(".tp"):
        raise ValueError(f"Unexpected extension: {args.infile}")

    if not args.outdir:
        raise ValueError("Missing `--outdir` argument to tempita.py")

    outdir_abs = os.path.join(os.getcwd(), args.outdir)
    outfile = os.path.join(
        outdir_abs, os.path.splitext(os.path.split(args.infile)[1])[0]
    )

    process_tempita(args.infile, outfile)


if __name__ == "__main__":
    main()
