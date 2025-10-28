#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import shutil
import sys
from pathlib import Path
import typing as T

def run(argsv: T.List[str]) -> int:
    commands: T.List[T.List[str]] = [[]]
    SEPARATOR = ';;;'

    # Generate CMD parameters
    parser = argparse.ArgumentParser(description='Wrapper for add_custom_command')
    parser.add_argument('-d', '--directory', type=str, metavar='D', required=True, help='Working directory to cwd to')
    parser.add_argument('-o', '--outputs', nargs='+', metavar='O', required=True, help='Expected output files')
    parser.add_argument('-O', '--original-outputs', nargs='*', metavar='O', default=[], help='Output files expected by CMake')
    parser.add_argument('commands', nargs=argparse.REMAINDER, help=f'A "{SEPARATOR}" separated list of commands')

    # Parse
    args = parser.parse_args(argsv)
    directory = Path(args.directory)

    dummy_target = None
    if len(args.outputs) == 1 and len(args.original_outputs) == 0:
        dummy_target = Path(args.outputs[0])
    elif len(args.outputs) != len(args.original_outputs):
        print('Length of output list and original output list differ')
        return 1

    for i in args.commands:
        if i == SEPARATOR:
            commands += [[]]
            continue

        i = i.replace('"', '')  # Remove leftover quotes
        commands[-1] += [i]

    # Execute
    for i in commands:
        # Skip empty lists
        if not i:
            continue

        cmd = []
        stdout = None
        stderr = None
        capture_file = ''

        for j in i:
            if j in {'>', '>>'}:
                stdout = subprocess.PIPE
                continue
            elif j in {'&>', '&>>'}:
                stdout = subprocess.PIPE
                stderr = subprocess.STDOUT
                continue

            if stdout is not None or stderr is not None:
                capture_file += j
            else:
                cmd += [j]

        try:
            directory.mkdir(parents=True, exist_ok=True)

            res = subprocess.run(cmd, stdout=stdout, stderr=stderr, cwd=str(directory), check=True)
            if capture_file:
                out_file = directory / capture_file
                out_file.write_bytes(res.stdout)
        except subprocess.CalledProcessError:
            return 1

    if dummy_target:
        dummy_target.touch()
        return 0

    # Copy outputs
    zipped_outputs = zip([Path(x) for x in args.outputs], [Path(x) for x in args.original_outputs])
    for expected, generated in zipped_outputs:
        do_copy = False
        if not expected.exists():
            if not generated.exists():
                print('Unable to find generated file. This can cause the build to fail:')
                print(generated)
                do_copy = False
            else:
                do_copy = True
        elif generated.exists():
            if generated.stat().st_mtime > expected.stat().st_mtime:
                do_copy = True

        if do_copy:
            if expected.exists():
                expected.unlink()
            shutil.copyfile(str(generated), str(expected))

    return 0

if __name__ == '__main__':
    sys.exit(run(sys.argv[1:]))
