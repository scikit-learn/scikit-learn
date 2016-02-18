#!/usr/bin/env python
""" cythonize

Cythonize pyx files into C files as needed.

Usage: cythonize [root_dir]

Default [root_dir] is 'sklearn'.

Checks pyx files to see if they have been changed relative to their
corresponding C files.  If they have, then runs cython on these files to
recreate the C files.

The script detects changes in the pyx/pxd files using checksums
[or hashes] stored in a database file

Simple script to invoke Cython on all .pyx
files; while waiting for a proper build system. Uses file hashes to
figure out if rebuild is needed.

It is called by ./setup.py sdist so that sdist package can be installed without
cython

Originally written by Dag Sverre Seljebotn, and adapted from statsmodel 0.6.1
(Modified BSD 3-clause)

We copied it for scikit-learn.

Note: this script does not check any of the dependent C libraries; it only
operates on the Cython .pyx files or their corresponding Cython header (.pxd)
files.
"""
# Author: Arthur Mensch <arthur.mensch@inria.fr>
# Author: Raghav R V <rvraghav93@gmail.com>
#
# License: BSD 3 clause

from __future__ import division, print_function, absolute_import

import os
import re
import sys
import hashlib
import subprocess

HASH_FILE = 'cythonize.dat'
DEFAULT_ROOT = 'sklearn'

# WindowsError is not defined on unix systems
try:
    WindowsError
except NameError:
    WindowsError = None


def cythonize(cython_file, gen_file):
    try:
        from Cython.Compiler.Version import version as cython_version
        from distutils.version import LooseVersion
        if LooseVersion(cython_version) < LooseVersion('0.21'):
            raise Exception('Building scikit-learn requires Cython >= 0.21')

    except ImportError:
        pass

    flags = ['--fast-fail']
    if gen_file.endswith('.cpp'):
        flags += ['--cplus']

    try:
        try:
            rc = subprocess.call(['cython'] +
                                 flags + ["-o", gen_file, cython_file])
            if rc != 0:
                raise Exception('Cythonizing %s failed' % cython_file)
        except OSError:
            # There are ways of installing Cython that don't result in a cython
            # executable on the path, see scipy issue gh-2397.
            rc = subprocess.call([sys.executable, '-c',
                                  'import sys; from Cython.Compiler.Main '
                                  'import setuptools_main as main;'
                                  ' sys.exit(main())'] + flags +
                                 ["-o", gen_file, cython_file])
            if rc != 0:
                raise Exception('Cythonizing %s failed' % cython_file)
    except OSError:
        raise OSError('Cython needs to be installed')


def load_hashes(filename):
    """Load the hashes dict from the hashfile"""
    # { filename : (sha1 of header if available or 'NA',
    #               sha1 of input,
    #               sha1 of output) }

    hashes = {}
    try:
        with open(filename, 'r') as cython_hash_file:
            for hash_record in cython_hash_file:
                (filename, header_hash,
                 cython_hash, gen_file_hash) = hash_record.split()
                hashes[filename] = (header_hash, cython_hash, gen_file_hash)
    except (KeyError, ValueError, AttributeError, IOError):
        hashes = {}
    return hashes


def save_hashes(hashes, filename):
    """Save the hashes dict to the hashfile"""
    with open(filename, 'w') as cython_hash_file:
        for key, value in hashes.items():
            cython_hash_file.write("%s %s %s %s\n"
                                   % (key, value[0], value[1], value[2]))


def sha1_of_file(filename):
    h = hashlib.sha1()
    with open(filename, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def clean_path(path):
    """Clean the path"""
    path = path.replace(os.sep, '/')
    if path.startswith('./'):
        path = path[2:]
    return path


def get_hash_tuple(header_path, cython_path, gen_file_path):
    """Get the hashes from the given files"""

    header_hash = (sha1_of_file(header_path)
                   if os.path.exists(header_path) else 'NA')
    from_hash = sha1_of_file(cython_path)
    to_hash = (sha1_of_file(gen_file_path)
               if os.path.exists(gen_file_path) else 'NA')

    return header_hash, from_hash, to_hash


def cythonize_if_unchanged(path, cython_file, gen_file, hashes):
    full_cython_path = os.path.join(path, cython_file)
    full_header_path = full_cython_path.replace('.pyx', '.pxd')
    full_gen_file_path = os.path.join(path, gen_file)

    current_hash = get_hash_tuple(full_header_path, full_cython_path,
                                  full_gen_file_path)

    if current_hash == hashes.get(clean_path(full_cython_path)):
        print('%s has not changed' % full_cython_path)
        return

    print('Processing %s' % full_cython_path)
    cythonize(full_cython_path, full_gen_file_path)

    # changed target file, recompute hash
    current_hash = get_hash_tuple(full_header_path, full_cython_path,
                                  full_gen_file_path)

    # Update the hashes dict with the new hash
    hashes[clean_path(full_cython_path)] = current_hash


def check_and_cythonize(root_dir):
    print(root_dir)
    hashes = load_hashes(HASH_FILE)

    for cur_dir, dirs, files in os.walk(root_dir):
        for filename in files:
            if filename.endswith('.pyx'):
                gen_file_ext = '.c'
                # Cython files with libcpp imports should be compiled to cpp
                with open(os.path.join(cur_dir, filename), 'rb') as f:
                    data = f.read()
                    m = re.search(b"libcpp", data, re.I | re.M)
                    if m:
                        gen_file_ext = ".cpp"
                cython_file = filename
                gen_file = filename.replace('.pyx', gen_file_ext)
                cythonize_if_unchanged(cur_dir, cython_file, gen_file, hashes)

                # Save hashes once per module. This prevents cythonizing prev.
                # files again when debugging broken code in a single file
                save_hashes(hashes, HASH_FILE)


def main(root_dir=DEFAULT_ROOT):
    check_and_cythonize(root_dir)


if __name__ == '__main__':
    try:
        root_dir_arg = sys.argv[1]
    except IndexError:
        root_dir_arg = DEFAULT_ROOT
    main(root_dir_arg)
