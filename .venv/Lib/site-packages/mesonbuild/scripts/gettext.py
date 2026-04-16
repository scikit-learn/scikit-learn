# SPDX-License-Identifier: Apache-2.0
# Copyright 2016 The Meson development team

from __future__ import annotations

import os
import argparse
import subprocess
import typing as T

parser = argparse.ArgumentParser()
parser.add_argument('command')
parser.add_argument('--pkgname', default='')
parser.add_argument('--datadirs', default='')
parser.add_argument('--langs', default='')
parser.add_argument('--localedir', default='')
parser.add_argument('--source-root', default='')
parser.add_argument('--subdir', default='')
parser.add_argument('--xgettext', default='xgettext')
parser.add_argument('--msgmerge', default='msgmerge')
parser.add_argument('--msginit', default='msginit')
parser.add_argument('--extra-args', default='')

def read_linguas(src_sub: str) -> T.List[str]:
    # Syntax of this file is documented here:
    # https://www.gnu.org/software/gettext/manual/html_node/po_002fLINGUAS.html
    linguas = os.path.join(src_sub, 'LINGUAS')
    try:
        langs = []
        with open(linguas, encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    langs += line.split()
        return langs
    except (FileNotFoundError, PermissionError):
        print(f'Could not find file LINGUAS in {src_sub}')
        return []

def run_potgen(src_sub: str, xgettext: str, pkgname: str, datadirs: str, args: T.List[str], source_root: str) -> int:
    listfile = os.path.join(src_sub, 'POTFILES.in')
    if not os.path.exists(listfile):
        listfile = os.path.join(src_sub, 'POTFILES')
        if not os.path.exists(listfile):
            print('Could not find file POTFILES in %s' % src_sub)
            return 1

    child_env = os.environ.copy()
    if datadirs:
        child_env['GETTEXTDATADIRS'] = datadirs

    ofile = os.path.join(src_sub, pkgname + '.pot')
    return subprocess.call([xgettext, '--package-name=' + pkgname, '-p', src_sub, '-f', listfile,
                            '-D', source_root, '-k_', '-o', ofile] + args,
                           env=child_env)

def update_po(src_sub: str, msgmerge: str, msginit: str, pkgname: str, langs: T.List[str]) -> int:
    potfile = os.path.join(src_sub, pkgname + '.pot')
    for l in langs:
        pofile = os.path.join(src_sub, l + '.po')
        if os.path.exists(pofile):
            subprocess.check_call([msgmerge, '-q', '-o', pofile, pofile, potfile])
        else:
            subprocess.check_call([msginit, '--input', potfile, '--output-file', pofile, '--locale', l, '--no-translator'])
    return 0

def run(args: T.List[str]) -> int:
    options = parser.parse_args(args)
    subcmd = options.command
    langs = options.langs.split('@@') if options.langs else None
    extra_args = options.extra_args.split('@@') if options.extra_args else []
    subdir = options.subdir
    src_sub = os.path.join(options.source_root, subdir)

    if not langs:
        langs = read_linguas(src_sub)

    if subcmd == 'pot':
        return run_potgen(src_sub, options.xgettext, options.pkgname, options.datadirs, extra_args, options.source_root)
    elif subcmd == 'update_po':
        if run_potgen(src_sub, options.xgettext, options.pkgname, options.datadirs, extra_args, options.source_root) != 0:
            return 1
        return update_po(src_sub, options.msgmerge, options.msginit, options.pkgname, langs)
    else:
        print('Unknown subcommand.')
        return 1
