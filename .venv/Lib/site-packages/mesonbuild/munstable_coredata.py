# SPDX-License-Identifier: Apache-2.0
# Copyright 2019 The Meson development team

from __future__ import annotations


from . import coredata as cdata
from .mesonlib import MachineChoice
from .options import OptionKey

import os.path
import pprint
import textwrap

# Note: when adding arguments, please also add them to the completion
# scripts in $MESONSRC/data/shell-completions/
def add_arguments(parser):
    parser.add_argument('--all', action='store_true', dest='all', default=False,
                        help='Show data not used by current backend.')

    parser.add_argument('builddir', nargs='?', default='.', help='The build directory')


def dump_compilers(compilers):
    for lang, compiler in compilers.items():
        print('  ' + lang + ':')
        print('      Id: ' + compiler.id)
        print('      Command: ' + ' '.join(compiler.exelist))
        if compiler.full_version:
            print('      Full version: ' + compiler.full_version)
        if compiler.version:
            print('      Detected version: ' + compiler.version)


def dump_guids(d):
    for name, value in d.items():
        print('  ' + name + ': ' + value)


def run(options):
    datadir = 'meson-private'
    if options.builddir is not None:
        datadir = os.path.join(options.builddir, datadir)
    if not os.path.isdir(datadir):
        print('Current directory is not a build dir. Please specify it or '
              'change the working directory to it.')
        return 1

    all_backends = options.all

    print('This is a dump of the internal unstable cache of meson. This is for debugging only.')
    print('Do NOT parse, this will change from version to version in incompatible ways')
    print('')

    coredata = cdata.load(options.builddir)
    backend = coredata.optstore.get_value_for(OptionKey('backend'))
    for k, v in sorted(coredata.__dict__.items()):
        if k in {'backend_options', 'base_options', 'builtins', 'compiler_options', 'user_options'}:
            # use `meson configure` to view these
            pass
        elif k in {'install_guid', 'test_guid', 'regen_guid'}:
            if all_backends or backend.startswith('vs'):
                print(k + ': ' + v)
        elif k == 'target_guids':
            if all_backends or backend.startswith('vs'):
                print(k + ':')
                dump_guids(v)
        elif k == 'lang_guids':
            if all_backends or backend.startswith('vs') or backend == 'xcode':
                print(k + ':')
                dump_guids(v)
        elif k == 'meson_command':
            if all_backends or backend.startswith('vs'):
                print('Meson command used in build file regeneration: ' + ' '.join(v))
        elif k == 'pkgconf_envvar':
            print('Last seen PKGCONFIG environment variable value: ' + v)
        elif k == 'version':
            print('Meson version: ' + v)
        elif k == 'cross_files':
            if v:
                print('Cross File: ' + ' '.join(v))
        elif k == 'config_files':
            if v:
                print('Native File: ' + ' '.join(v))
        elif k == 'compilers':
            for for_machine in MachineChoice:
                print('Cached {} machine compilers:'.format(
                    for_machine.get_lower_case_name()))
                dump_compilers(v[for_machine])
        elif k == 'deps':
            def print_dep(dep_key, dep):
                print('  ' + dep_key[0][1] + ": ")
                print('      compile args: ' + repr(dep.get_compile_args()))
                print('      link args: ' + repr(dep.get_link_args()))
                if dep.get_sources():
                    print('      sources: ' + repr(dep.get_sources()))
                print('      version: ' + repr(dep.get_version()))

            for for_machine in iter(MachineChoice):
                items_list = sorted(v[for_machine].items())
                if items_list:
                    print(f'Cached dependencies for {for_machine.get_lower_case_name()} machine')
                    for dep_key, deps in items_list:
                        for dep in deps:
                            print_dep(dep_key, dep)
        else:
            print(k + ':')
            print(textwrap.indent(pprint.pformat(v), '  '))
