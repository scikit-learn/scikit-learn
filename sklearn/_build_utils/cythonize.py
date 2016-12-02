#!/usr/bin/env python
""" cythonize

Cythonize pyx files into C files as needed.

Usage: cythonize [root_dir]

Default [root_dir] is 'sklearn'.

Checks pyx files to see if they have been changed relative to their
corresponding C files.  If they have, then runs cython on these files to
recreate the C files.

The script detects changes in the pyx files using checksums [or hashes] stored
in a database file

Simple script to invoke Cython on all .pyx
files; while waiting for a proper build system. Uses file hashes to
figure out if rebuild is needed.

It is called by ./setup.py sdist so that sdist package can be installed without
cython

Originally written by Dag Sverre Seljebotn, and adapted from statsmodel 0.6.1
(Modified BSD 3-clause)

We copied it for scikit-learn.

Note: this script does not check any of the dependent C libraries; it only
operates on the Cython .pyx files.
"""
# author: Arthur Mensch
# license: BSD

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


#
# Rules
#
def process_pyx(fromfile, tofile):
    try:
        from Cython.Compiler.Version import version as cython_version
        from distutils.version import LooseVersion
        if LooseVersion(cython_version) < LooseVersion('0.21'):
            raise Exception('Building scikit-learn requires Cython >= 0.21')

    except ImportError:
        pass

    flags = ['--fast-fail']
    if tofile.endswith('.cpp'):
        flags += ['--cplus']

    try:
        try:
            r = subprocess.call(['cython'] + flags + ["-o", tofile, fromfile])
            if r != 0:
                raise Exception('Cython failed')
        except OSError:
            # There are ways of installing Cython that don't result in a cython
            # executable on the path, see gh-2397.
            r = subprocess.call([sys.executable, '-c',
                                 'import sys; from Cython.Compiler.Main '
                                 'import setuptools_main as main;'
                                 ' sys.exit(main())'] + flags +
                                ["-o", tofile, fromfile])
            if r != 0:
                raise Exception('Cython failed')
    except OSError:
        raise OSError('Cython needs to be installed')


rules = {
    '.pyx': process_pyx,
}


#
# Hash db
#
def load_hashes(filename):
    # Return { filename : (sha1 of input, sha1 of output) }
    if os.path.isfile(filename):
        hashes = {}
        with open(filename, 'r') as f:
            for line in f:
                filename, inhash, outhash = line.split()
                hashes[filename] = (inhash, outhash)
    else:
        hashes = {}
    return hashes


def save_hashes(hash_db, filename):
    with open(filename, 'w') as f:
        for key, value in hash_db.items():
            f.write("%s %s %s\n" % (key, value[0], value[1]))


def sha1_of_file(filename):
    h = hashlib.sha1()
    with open(filename, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


#
# Main program
#
def normpath(path):
    path = path.replace(os.sep, '/')
    if path.startswith('./'):
        path = path[2:]
    return path


def get_hash(frompath, topath):
    from_hash = sha1_of_file(frompath)
    to_hash = sha1_of_file(topath) if os.path.exists(topath) else None
    return from_hash, to_hash


def process(path, fromfile, tofile, processor_function, hash_db):
    fullfrompath = os.path.join(path, fromfile)
    fulltopath = os.path.join(path, tofile)
    current_hash = get_hash(fullfrompath, fulltopath)
    if current_hash == hash_db.get(normpath(fullfrompath)):
        print('%s has not changed' % fullfrompath)
        return

    print('Processing %s' % fullfrompath)
    processor_function(fullfrompath, fulltopath)
    # changed target file, recompute hash
    current_hash = get_hash(fullfrompath, fulltopath)
    # store hash in db
    hash_db[normpath(fullfrompath)] = current_hash


def find_process_files(root_dir):
    print(root_dir)
    hash_db = load_hashes(HASH_FILE)
    for cur_dir, dirs, files in os.walk(root_dir):
        for filename in files:
            for fromext, function in rules.items():
                if filename.endswith(fromext):
                    toext = ".c"
                    with open(os.path.join(cur_dir, filename), 'rb') as f:
                        data = f.read()
                        m = re.search(b"libcpp", data, re.I | re.M)
                        if m:
                            toext = ".cpp"
                    fromfile = filename
                    tofile = filename[:-len(fromext)] + toext
                    process(cur_dir, fromfile, tofile, function, hash_db)
    save_hashes(hash_db, HASH_FILE)


def main(root_dir=DEFAULT_ROOT):
    find_process_files(root_dir)


if __name__ == '__main__':
    try:
        root_dir_arg = sys.argv[1]
    except IndexError:
        root_dir_arg = DEFAULT_ROOT
    main(root_dir_arg)
