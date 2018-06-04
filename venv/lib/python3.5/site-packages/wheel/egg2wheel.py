#!/usr/bin/env python
import distutils.dist
import os.path
import re
import shutil
import sys
import tempfile
import zipfile
from argparse import ArgumentParser
from distutils.archive_util import make_archive
from glob import iglob

import wheel.bdist_wheel
from wheel.tool import WheelError
from wheel.wininst2wheel import _bdist_wheel_tag

egg_info_re = re.compile(r'''
    (?P<name>.+?)-(?P<ver>.+?)
    (-(?P<pyver>py\d\.\d)
     (-(?P<arch>.+?))?
    )?.egg$''', re.VERBOSE)


def egg2wheel(egg_path, dest_dir):
    filename = os.path.basename(egg_path)
    match = egg_info_re.match(filename)
    if not match:
        raise WheelError('Invalid egg file name: {}'.format(filename))

    egg_info = match.groupdict()
    dir = tempfile.mkdtemp(suffix="_e2w")
    if os.path.isfile(egg_path):
        # assume we have a bdist_egg otherwise
        egg = zipfile.ZipFile(egg_path)
        egg.extractall(dir)
    else:
        # support buildout-style installed eggs directories
        for pth in os.listdir(egg_path):
            src = os.path.join(egg_path, pth)
            if os.path.isfile(src):
                shutil.copy2(src, dir)
            else:
                shutil.copytree(src, os.path.join(dir, pth))

    pyver = egg_info['pyver']
    if pyver:
        pyver = pyver.replace('.', '')

    arch = (egg_info['arch'] or 'any').replace('.', '_').replace('-', '_')

    # assume all binary eggs are for CPython
    abi = 'cp' + pyver[2:] if arch != 'any' else 'none'

    root_is_purelib = egg_info['arch'] is None
    if root_is_purelib:
        bw = wheel.bdist_wheel.bdist_wheel(distutils.dist.Distribution())
    else:
        bw = _bdist_wheel_tag(distutils.dist.Distribution())

    bw.root_is_pure = root_is_purelib
    bw.python_tag = pyver
    bw.plat_name_supplied = True
    bw.plat_name = egg_info['arch'] or 'any'
    if not root_is_purelib:
        bw.full_tag_supplied = True
        bw.full_tag = (pyver, abi, arch)

    dist_info_dir = os.path.join(dir, '{name}-{ver}.dist-info'.format(**egg_info))
    bw.egg2dist(os.path.join(dir, 'EGG-INFO'), dist_info_dir)
    bw.write_wheelfile(dist_info_dir, generator='egg2wheel')
    bw.write_record(dir, dist_info_dir)
    wheel_name = '{name}-{ver}-{pyver}-{}-{}'.format(abi, arch, **egg_info)
    filename = make_archive(os.path.join(dest_dir, wheel_name), 'zip', root_dir=dir)
    os.rename(filename, filename[:-3] + 'whl')
    shutil.rmtree(dir)


def main():
    parser = ArgumentParser()
    parser.add_argument('eggs', nargs='*', help="Eggs to convert")
    parser.add_argument('--dest-dir', '-d', default=os.path.curdir,
                        help="Directory to store wheels (default %(default)s)")
    parser.add_argument('--verbose', '-v', action='store_true')
    args = parser.parse_args()
    for pat in args.eggs:
        for egg in iglob(pat):
            if args.verbose:
                print("{}... ".format(egg))
                sys.stdout.flush()

            egg2wheel(egg, args.dest_dir)
            if args.verbose:
                print("OK")


if __name__ == "__main__":
    main()
