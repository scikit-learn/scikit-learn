# SPDX-License-Identifier: Apache-2.0
# Copyright 2014-2016 The Meson development team
# Copyright © 2023-2024 Intel Corporation

from __future__ import annotations

import itertools
import hashlib
import shutil
import os
import textwrap
import typing as T
import collections

from . import build
from . import coredata
from . import options
from . import environment
from . import mesonlib
from . import mintro
from . import mlog
from .ast import AstIDGenerator, IntrospectionInterpreter
from .mesonlib import MachineChoice
from .options import OptionKey
from .optinterpreter import OptionInterpreter

if T.TYPE_CHECKING:
    from typing_extensions import Protocol
    from typing import Any
    from .options import UserOption
    import argparse

    class CMDOptions(coredata.SharedCMDOptions, Protocol):

        builddir: str
        clearcache: bool
        pager: bool

    # cannot be TV_Loggable, because non-ansidecorators do direct string concat
    LOGLINE = T.Union[str, mlog.AnsiDecorator]

# Note: when adding arguments, please also add them to the completion
# scripts in $MESONSRC/data/shell-completions/
def add_arguments(parser: 'argparse.ArgumentParser') -> None:
    coredata.register_builtin_arguments(parser)
    parser.add_argument('builddir', nargs='?', default='.')
    parser.add_argument('--clearcache', action='store_true', default=False,
                        help='Clear cached state (e.g. found dependencies)')
    parser.add_argument('--no-pager', action='store_false', dest='pager',
                        help='Do not redirect output to a pager')

def stringify(val: T.Any) -> str:
    if isinstance(val, bool):
        return str(val).lower()
    elif isinstance(val, list):
        s = ', '.join(stringify(i) for i in val)
        return f'[{s}]'
    elif val is None:
        return ''
    else:
        return str(val)


class ConfException(mesonlib.MesonException):
    pass


class Conf:
    def __init__(self, build_dir: str):
        self.build_dir = os.path.abspath(os.path.realpath(build_dir))
        if 'meson.build' in [os.path.basename(self.build_dir), self.build_dir]:
            self.build_dir = os.path.dirname(self.build_dir)
        self.build = None
        self.max_choices_line_length = 60
        self.name_col: T.List[LOGLINE] = []
        self.value_col: T.List[LOGLINE] = []
        self.choices_col: T.List[LOGLINE] = []
        self.descr_col: T.List[LOGLINE] = []
        self.all_subprojects: T.Set[str] = set()

        if os.path.isdir(os.path.join(self.build_dir, 'meson-private')):
            self.build = build.load(self.build_dir)
            self.source_dir = self.build.environment.get_source_dir()
            self.coredata = self.build.environment.coredata
            self.default_values_only = False

            # if the option file has been updated, reload it
            # This cannot handle options for a new subproject that has not yet
            # been configured.
            for sub, conf_options in self.coredata.options_files.items():
                if conf_options is not None and os.path.exists(conf_options[0]):
                    opfile = conf_options[0]
                    with open(opfile, 'rb') as f:
                        ophash = hashlib.sha1(f.read()).hexdigest()
                        if ophash != conf_options[1]:
                            oi = OptionInterpreter(self.coredata.optstore, sub)
                            oi.process(opfile)
                            self.coredata.update_project_options(oi.options, sub)
                            self.coredata.options_files[sub] = (opfile, ophash)
                else:
                    opfile = os.path.join(self.source_dir, 'meson.options')
                    if not os.path.exists(opfile):
                        opfile = os.path.join(self.source_dir, 'meson_options.txt')
                    if os.path.exists(opfile):
                        oi = OptionInterpreter(self.coredata.optstore, sub)
                        oi.process(opfile)
                        self.coredata.update_project_options(oi.options, sub)
                        with open(opfile, 'rb') as f:
                            ophash = hashlib.sha1(f.read()).hexdigest()
                        self.coredata.options_files[sub] = (opfile, ophash)
                    else:
                        self.coredata.update_project_options({}, sub)
        elif os.path.isfile(os.path.join(self.build_dir, environment.build_filename)):
            # Make sure that log entries in other parts of meson don't interfere with the JSON output
            with mlog.no_logging():
                self.source_dir = os.path.abspath(os.path.realpath(self.build_dir))
                intr = IntrospectionInterpreter(self.source_dir, '', 'ninja', visitors = [AstIDGenerator()])
                intr.analyze()
            self.coredata = intr.coredata
            self.default_values_only = True
        else:
            raise ConfException(f'Directory {build_dir} is neither a Meson build directory nor a project source directory.')

    def clear_cache(self) -> None:
        self.coredata.clear_cache()

    def set_options(self, options: T.Dict[OptionKey, str]) -> bool:
        return self.coredata.set_options(options)

    def save(self) -> None:
        # Do nothing when using introspection
        if self.default_values_only:
            return
        coredata.save(self.coredata, self.build_dir)
        # We don't write the build file because any changes to it
        # are erased when Meson is executed the next time, i.e. when
        # Ninja is run.

    def print_aligned(self) -> None:
        """Do the actual printing.

        This prints the generated output in an aligned, pretty form. it aims
        for a total width of 160 characters, but will use whatever the tty
        reports its value to be. Though this is much wider than the standard
        80 characters of terminals, and even than the newer 120, compressing
        it to those lengths makes the output hard to read.

        Each column will have a specific width, and will be line wrapped.
        """
        total_width = shutil.get_terminal_size(fallback=(160, 0))[0]
        _col = max(total_width // 5, 20)
        last_column = total_width - (3 * _col) - 3
        four_column = (_col, _col, _col, last_column if last_column > 1 else _col)

        for line in zip(self.name_col, self.value_col, self.choices_col, self.descr_col):
            if not any(line):
                mlog.log('')
                continue

            # This is a header, like `Subproject foo:`,
            # We just want to print that and get on with it
            if line[0] and not any(line[1:]):
                mlog.log(line[0])
                continue

            def wrap_text(text: LOGLINE, width: int) -> mlog.TV_LoggableList:
                raw = text.text if isinstance(text, mlog.AnsiDecorator) else text
                indent = ' ' if raw.startswith('[') else ''
                wrapped_ = textwrap.wrap(raw, width, subsequent_indent=indent)
                # We cast this because https://github.com/python/mypy/issues/1965
                # mlog.TV_LoggableList does not provide __len__ for stringprotocol
                if isinstance(text, mlog.AnsiDecorator):
                    wrapped = T.cast('T.List[LOGLINE]', [mlog.AnsiDecorator(i, text.code) for i in wrapped_])
                else:
                    wrapped = T.cast('T.List[LOGLINE]', wrapped_)
                # Add padding here to get even rows, as `textwrap.wrap()` will
                # only shorten, not lengthen each item
                return [str(i) + ' ' * (width - len(i)) for i in wrapped]

            # wrap will take a long string, and create a list of strings no
            # longer than the size given. Then that list can be zipped into, to
            # print each line of the output, such the that columns are printed
            # to the right width, row by row.
            name = wrap_text(line[0], four_column[0])
            val = wrap_text(line[1], four_column[1])
            choice = wrap_text(line[2], four_column[2])
            desc = wrap_text(line[3], four_column[3])
            for l in itertools.zip_longest(name, val, choice, desc, fillvalue=''):
                items = [l[i] if l[i] else ' ' * four_column[i] for i in range(4)]
                mlog.log(*items)

    def split_options_per_subproject(self, options: 'T.Union[dict[OptionKey, UserOption[Any]], coredata.KeyedOptionDictType]') -> T.Dict[str, 'coredata.MutableKeyedOptionDictType']:
        result: T.Dict[str, 'coredata.MutableKeyedOptionDictType'] = {}
        for k, o in options.items():
            if k.subproject:
                self.all_subprojects.add(k.subproject)
            result.setdefault(k.subproject, {})[k] = o
        return result

    def _add_line(self, name: LOGLINE, value: LOGLINE, choices: LOGLINE, descr: LOGLINE) -> None:
        if isinstance(name, mlog.AnsiDecorator):
            name.text = ' ' * self.print_margin + name.text
        else:
            name = ' ' * self.print_margin + name
        self.name_col.append(name)
        self.value_col.append(value)
        self.choices_col.append(choices)
        self.descr_col.append(descr)

    def add_option(self, name: str, descr: str, value: T.Any, choices: T.Any) -> None:
        value = stringify(value)
        choices = stringify(choices)
        self._add_line(mlog.green(name), mlog.yellow(value), mlog.blue(choices), descr)

    def add_title(self, title: str) -> None:
        newtitle = mlog.cyan(title)
        descr = mlog.cyan('Description')
        value = mlog.cyan('Default Value' if self.default_values_only else 'Current Value')
        choices = mlog.cyan('Possible Values')
        self._add_line('', '', '', '')
        self._add_line(newtitle, value, choices, descr)
        self._add_line('-' * len(newtitle), '-' * len(value), '-' * len(choices), '-' * len(descr))

    def add_section(self, section: str) -> None:
        self.print_margin = 0
        self._add_line('', '', '', '')
        self._add_line(mlog.normal_yellow(section + ':'), '', '', '')
        self.print_margin = 2

    def print_options(self, title: str, opts: 'T.Union[dict[OptionKey, UserOption[Any]], coredata.KeyedOptionDictType]') -> None:
        if not opts:
            return
        if title:
            self.add_title(title)
        auto = T.cast('options.UserFeatureOption', self.coredata.optstore.get_value_object('auto_features'))
        for k, o in sorted(opts.items()):
            printable_value = o.printable_value()
            root = k.as_root()
            if o.yielding and k.subproject and root in self.coredata.optstore:
                printable_value = '<inherited from main project>'
            if isinstance(o, options.UserFeatureOption) and o.is_auto():
                printable_value = auto.printable_value()
            self.add_option(str(root), o.description, printable_value, o.choices)

    def print_conf(self, pager: bool) -> None:
        if pager:
            mlog.start_pager()

        def print_default_values_warning() -> None:
            mlog.warning('The source directory instead of the build directory was specified.')
            mlog.warning('Only the default values for the project are printed.')

        if self.default_values_only:
            print_default_values_warning()
            mlog.log('')

        mlog.log('Core properties:')
        mlog.log('  Source dir', self.source_dir)
        if not self.default_values_only:
            mlog.log('  Build dir ', self.build_dir)

        dir_option_names = set(options.BUILTIN_DIR_OPTIONS)
        test_option_names = {OptionKey('errorlogs'),
                             OptionKey('stdsplit')}

        dir_options: 'coredata.MutableKeyedOptionDictType' = {}
        test_options: 'coredata.MutableKeyedOptionDictType' = {}
        core_options: 'coredata.MutableKeyedOptionDictType' = {}
        module_options: T.Dict[str, 'coredata.MutableKeyedOptionDictType'] = collections.defaultdict(dict)
        for k, v in self.coredata.optstore.items():
            if k in dir_option_names:
                dir_options[k] = v
            elif k in test_option_names:
                test_options[k] = v
            elif k.has_module_prefix():
                # Ignore module options if we did not use that module during
                # configuration.
                modname = k.get_module_prefix()
                if self.build and modname not in self.build.modules:
                    continue
                module_options[modname][k] = v
            elif self.coredata.optstore.is_builtin_option(k):
                core_options[k] = v

        host_core_options = self.split_options_per_subproject({k: v for k, v in core_options.items() if k.machine is MachineChoice.HOST})
        build_core_options = self.split_options_per_subproject({k: v for k, v in core_options.items() if k.machine is MachineChoice.BUILD})
        host_compiler_options = self.split_options_per_subproject({k: v for k, v in self.coredata.optstore.items() if self.coredata.optstore.is_compiler_option(k) and k.machine is MachineChoice.HOST})
        build_compiler_options = self.split_options_per_subproject({k: v for k, v in self.coredata.optstore.items() if self.coredata.optstore.is_compiler_option(k) and k.machine is MachineChoice.BUILD})
        project_options = self.split_options_per_subproject({k: v for k, v in self.coredata.optstore.items() if self.coredata.optstore.is_project_option(k)})
        show_build_options = self.default_values_only or self.build.environment.is_cross_build()

        self.add_section('Main project options')
        self.print_options('Core options', host_core_options[''])
        if show_build_options:
            self.print_options('', build_core_options[''])
        self.print_options('Backend options', {k: v for k, v in self.coredata.optstore.items() if self.coredata.optstore.is_backend_option(k)})
        self.print_options('Base options', {k: v for k, v in self.coredata.optstore.items() if self.coredata.optstore.is_base_option(k)})
        self.print_options('Compiler options', host_compiler_options.get('', {}))
        if show_build_options:
            self.print_options('', build_compiler_options.get('', {}))
        for mod, mod_options in module_options.items():
            self.print_options(f'{mod} module options', mod_options)
        self.print_options('Directories', dir_options)
        self.print_options('Testing options', test_options)
        self.print_options('Project options', project_options.get('', {}))
        for subproject in sorted(self.all_subprojects):
            if subproject == '':
                continue
            self.add_section('Subproject ' + subproject)
            if subproject in host_core_options:
                self.print_options('Core options', host_core_options[subproject])
            if subproject in build_core_options and show_build_options:
                self.print_options('', build_core_options[subproject])
            if subproject in host_compiler_options:
                self.print_options('Compiler options', host_compiler_options[subproject])
            if subproject in build_compiler_options and show_build_options:
                self.print_options('', build_compiler_options[subproject])
            if subproject in project_options:
                self.print_options('Project options', project_options[subproject])
        self.print_aligned()

        # Print the warning twice so that the user shouldn't be able to miss it
        if self.default_values_only:
            mlog.log('')
            print_default_values_warning()

        self.print_nondefault_buildtype_options()

    def print_nondefault_buildtype_options(self) -> None:
        mismatching = self.coredata.get_nondefault_buildtype_args()
        if not mismatching:
            return
        mlog.log("\nThe following option(s) have a different value than the build type default\n")
        mlog.log('               current   default')
        for m in mismatching:
            mlog.log(f'{m[0]:21}{m[1]:10}{m[2]:10}')

def run_impl(options: CMDOptions, builddir: str) -> int:
    print_only = not options.cmd_line_options and not options.clearcache
    c = None
    try:
        c = Conf(builddir)
        if c.default_values_only and not print_only:
            raise mesonlib.MesonException('No valid build directory found, cannot modify options.')
        if c.default_values_only or print_only:
            c.print_conf(options.pager)
            return 0

        save = False
        if options.cmd_line_options:
            save = c.set_options(options.cmd_line_options)
            coredata.update_cmd_line_file(builddir, options)
        if options.clearcache:
            c.clear_cache()
            save = True
        if save:
            c.save()
            mintro.update_build_options(c.coredata, c.build.environment.info_dir)
            mintro.write_meson_info_file(c.build, [])
    except ConfException as e:
        mlog.log('Meson configurator encountered an error:')
        if c is not None and c.build is not None:
            mintro.write_meson_info_file(c.build, [e])
        raise e
    except BrokenPipeError:
        # Pager quit before we wrote everything.
        pass
    return 0

def run(options: CMDOptions) -> int:
    coredata.parse_cmd_line_options(options)
    builddir = os.path.abspath(os.path.realpath(options.builddir))
    return run_impl(options, builddir)
