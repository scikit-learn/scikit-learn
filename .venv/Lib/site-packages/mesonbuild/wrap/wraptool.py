# SPDX-License-Identifier: Apache-2.0
# Copyright 2015-2016 The Meson development team

from __future__ import annotations

import sys, os
import configparser
import shutil
import typing as T

from glob import glob
from .wrap import (open_wrapdburl, read_and_decompress, WrapException, get_releases,
                   get_releases_data, parse_patch_url)
from .. import mesonlib, msubprojects

if T.TYPE_CHECKING:
    import argparse

# Note: when adding arguments, please also add them to the completion
# scripts in $MESONSRC/data/shell-completions/
def add_arguments(parser: 'argparse.ArgumentParser') -> None:
    subparsers = parser.add_subparsers(title='Commands', dest='command')
    subparsers.required = True

    p = subparsers.add_parser('list', help='show all available projects')
    p.add_argument('--allow-insecure', default=False, action='store_true',
                   help='Allow insecure server connections.')
    p.set_defaults(wrap_func=list_projects)

    p = subparsers.add_parser('search', help='search the db by name')
    p.add_argument('--allow-insecure', default=False, action='store_true',
                   help='Allow insecure server connections.')
    p.add_argument('name')
    p.set_defaults(wrap_func=search)

    p = subparsers.add_parser('install', help='install the specified project')
    p.add_argument('--allow-insecure', default=False, action='store_true',
                   help='Allow insecure server connections.')
    p.add_argument('name')
    p.set_defaults(wrap_func=install)

    p = msubprojects.add_wrap_update_parser(subparsers)
    p.set_defaults(wrap_func=msubprojects.run)

    p = subparsers.add_parser('info', help='show available versions of a project')
    p.add_argument('--allow-insecure', default=False, action='store_true',
                   help='Allow insecure server connections.')
    p.add_argument('name')
    p.set_defaults(wrap_func=info)

    p = subparsers.add_parser('status', help='show installed and available versions of your projects')
    p.add_argument('--allow-insecure', default=False, action='store_true',
                   help='Allow insecure server connections.')
    p.set_defaults(wrap_func=status)

    p = subparsers.add_parser('promote', help='bring a subsubproject up to the master project')
    p.add_argument('project_path')
    p.set_defaults(wrap_func=promote)

    p = subparsers.add_parser('update-db', help='Update list of projects available in WrapDB (Since 0.61.0)')
    p.add_argument('--allow-insecure', default=False, action='store_true',
                   help='Allow insecure server connections.')
    p.set_defaults(wrap_func=update_db)

def list_projects(options: 'argparse.Namespace') -> None:
    releases = get_releases(options.allow_insecure)
    for p in releases.keys():
        print(p)

def search(options: 'argparse.Namespace') -> None:
    name = options.name
    releases = get_releases(options.allow_insecure)
    for p, info in releases.items():
        if p.find(name) != -1:
            print(p)
        else:
            for dep in info.get('dependency_names', []):
                if dep.find(name) != -1:
                    print(f'Dependency {dep} found in wrap {p}')

def get_latest_version(name: str, allow_insecure: bool) -> T.Tuple[str, str]:
    releases = get_releases(allow_insecure)
    info = releases.get(name)
    if not info:
        raise WrapException(f'Wrap {name} not found in wrapdb')
    latest_version = info['versions'][0]
    version, revision = latest_version.rsplit('-', 1)
    return version, revision

def install(options: 'argparse.Namespace') -> None:
    name = options.name
    subproject_dir_name = mesonlib.get_subproject_dir()
    if subproject_dir_name is None or not os.path.isdir(subproject_dir_name):
        raise SystemExit('Subprojects dir not found. Run this script in your source root directory.')
    if os.path.isdir(os.path.join(subproject_dir_name, name)):
        raise SystemExit('Subproject directory for this project already exists.')
    wrapfile = os.path.join(subproject_dir_name, name + '.wrap')
    if os.path.exists(wrapfile):
        raise SystemExit('Wrap file already exists.')
    (version, revision) = get_latest_version(name, options.allow_insecure)
    url = open_wrapdburl(f'https://wrapdb.mesonbuild.com/v2/{name}_{version}-{revision}/{name}.wrap', options.allow_insecure, True, True)
    with open(wrapfile, 'wb') as f:
        f.write(read_and_decompress(url))
    print(f'Installed {name} version {version} revision {revision}')

def get_current_version(wrapfile: str) -> T.Tuple[str, str, str, str, T.Optional[str]]:
    cp = configparser.ConfigParser(interpolation=None)
    cp.read(wrapfile)
    try:
        wrap_data = cp['wrap-file']
    except KeyError:
        raise WrapException('Not a wrap-file, cannot have come from the wrapdb')
    try:
        patch_url = wrap_data['patch_url']
    except KeyError:
        # We assume a wrap without a patch_url is probably just an pointer to upstream's
        # build files. The version should be in the tarball filename, even if it isn't
        # purely guaranteed. The wrapdb revision should be 1 because it just needs uploading once.
        branch = mesonlib.search_version(wrap_data['source_filename'])
        revision, patch_filename = '1', None
    else:
        branch, revision = parse_patch_url(patch_url)
        patch_filename = wrap_data['patch_filename']
    return branch, revision, wrap_data['directory'], wrap_data['source_filename'], patch_filename

def info(options: 'argparse.Namespace') -> None:
    name = options.name
    releases = get_releases(options.allow_insecure)
    info = releases.get(name)
    if not info:
        raise WrapException(f'Wrap {name} not found in wrapdb')
    print(f'Available versions of {name}:')
    for v in info['versions']:
        print(' ', v)

def do_promotion(from_path: str, spdir_name: str) -> None:
    if os.path.isfile(from_path):
        assert from_path.endswith('.wrap')
        shutil.copy(from_path, spdir_name)
    elif os.path.isdir(from_path):
        sproj_name = os.path.basename(from_path)
        outputdir = os.path.join(spdir_name, sproj_name)
        if os.path.exists(outputdir):
            raise SystemExit(f'Output dir {outputdir} already exists. Will not overwrite.')

        subpdir = mesonlib.get_subproject_dir()
        if subpdir is not None:
            ignore = shutil.ignore_patterns(subpdir)
        else:
            ignore = None

        shutil.copytree(from_path, outputdir, ignore=ignore)

def promote(options: 'argparse.Namespace') -> None:
    argument = options.project_path
    spdir_name = mesonlib.get_subproject_dir()
    if spdir_name is None:
        raise SystemExit('Subproject dir not found. Run this script in your source root directory.')
    sprojs = mesonlib.detect_subprojects(spdir_name)

    # check if the argument is a full path to a subproject directory or wrap file
    system_native_path_argument = argument.replace('/', os.sep)
    for matches in sprojs.values():
        if system_native_path_argument in matches:
            do_promotion(system_native_path_argument, spdir_name)
            return

    # otherwise the argument is just a subproject basename which must be unambiguous
    if argument not in sprojs:
        raise SystemExit(f'Subproject {argument} not found in directory tree.')
    matches = sprojs[argument]
    if len(matches) > 1:
        print(f'There is more than one version of {argument} in tree. Please specify which one to promote:\n', file=sys.stderr)
        for s in matches:
            print(s, file=sys.stderr)
        raise SystemExit(1)
    do_promotion(matches[0], spdir_name)

def status(options: 'argparse.Namespace') -> None:
    print('Subproject status')
    subdir = mesonlib.get_subproject_dir()
    assert subdir is not None, "This should only happen in a non-native subproject"
    for w in glob(f'{subdir}/*.wrap'):
        name = os.path.basename(w)[:-5]
        try:
            (latest_branch, latest_revision) = get_latest_version(name, options.allow_insecure)
        except Exception:
            print('', name, 'not available in wrapdb.', file=sys.stderr)
            continue
        try:
            (current_branch, current_revision, _, _, _) = get_current_version(w)
        except Exception:
            print('', name, 'Wrap file not from wrapdb.', file=sys.stderr)
            continue
        if current_branch == latest_branch and current_revision == latest_revision:
            print('', name, f'up to date. Branch {current_branch}, revision {current_revision}.')
        else:
            print('', name, f'not up to date. Have {current_branch} {current_revision}, but {latest_branch} {latest_revision} is available.')

def update_db(options: 'argparse.Namespace') -> None:
    data = get_releases_data(options.allow_insecure)
    subproject_dir_name = mesonlib.get_subproject_dir()
    if subproject_dir_name is None:
        raise SystemExit('Subproject dir not found. Run this script in your source root directory.')

    os.makedirs(subproject_dir_name, exist_ok=True)
    with open(os.path.join(subproject_dir_name, 'wrapdb.json'), 'wb') as f:
        f.write(data)

def run(options: 'argparse.Namespace') -> int:
    options.wrap_func(options)
    return 0
