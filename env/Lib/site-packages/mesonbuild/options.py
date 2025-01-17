# SPDX-License-Identifier: Apache-2.0
# Copyright 2013-2024 Contributors to the The Meson project
# Copyright © 2019-2024 Intel Corporation

from __future__ import annotations
from collections import OrderedDict
from itertools import chain
from functools import total_ordering
import argparse
import typing as T

from .mesonlib import (
    HoldableObject,
    default_prefix,
    default_datadir,
    default_includedir,
    default_infodir,
    default_libdir,
    default_libexecdir,
    default_localedir,
    default_mandir,
    default_sbindir,
    default_sysconfdir,
    MesonException,
    listify_array_value,
    MachineChoice,
)
from . import mlog

if T.TYPE_CHECKING:
    from typing_extensions import TypedDict

    class ArgparseKWs(TypedDict, total=False):

        action: str
        dest: str
        default: str
        choices: T.List

DEFAULT_YIELDING = False

# Can't bind this near the class method it seems, sadly.
_T = T.TypeVar('_T')

backendlist = ['ninja', 'vs', 'vs2010', 'vs2012', 'vs2013', 'vs2015', 'vs2017', 'vs2019', 'vs2022', 'xcode', 'none']
genvslitelist = ['vs2022']
buildtypelist = ['plain', 'debug', 'debugoptimized', 'release', 'minsize', 'custom']


# This is copied from coredata. There is no way to share this, because this
# is used in the OptionKey constructor, and the coredata lists are
# OptionKeys...
_BUILTIN_NAMES = {
    'prefix',
    'bindir',
    'datadir',
    'includedir',
    'infodir',
    'libdir',
    'licensedir',
    'libexecdir',
    'localedir',
    'localstatedir',
    'mandir',
    'sbindir',
    'sharedstatedir',
    'sysconfdir',
    'auto_features',
    'backend',
    'buildtype',
    'debug',
    'default_library',
    'default_both_libraries',
    'errorlogs',
    'genvslite',
    'install_umask',
    'layout',
    'optimization',
    'prefer_static',
    'stdsplit',
    'strip',
    'unity',
    'unity_size',
    'warning_level',
    'werror',
    'wrap_mode',
    'force_fallback_for',
    'pkg_config_path',
    'cmake_prefix_path',
    'vsenv',
}

@total_ordering
class OptionKey:

    """Represents an option key in the various option dictionaries.

    This provides a flexible, powerful way to map option names from their
    external form (things like subproject:build.option) to something that
    internally easier to reason about and produce.
    """

    __slots__ = ['name', 'subproject', 'machine', '_hash']

    name: str
    subproject: str
    machine: MachineChoice
    _hash: int

    def __init__(self, name: str, subproject: str = '',
                 machine: MachineChoice = MachineChoice.HOST):
        # the _type option to the constructor is kinda private. We want to be
        # able to save the state and avoid the lookup function when
        # pickling/unpickling, but we need to be able to calculate it when
        # constructing a new OptionKey
        object.__setattr__(self, 'name', name)
        object.__setattr__(self, 'subproject', subproject)
        object.__setattr__(self, 'machine', machine)
        object.__setattr__(self, '_hash', hash((name, subproject, machine)))

    def __setattr__(self, key: str, value: T.Any) -> None:
        raise AttributeError('OptionKey instances do not support mutation.')

    def __getstate__(self) -> T.Dict[str, T.Any]:
        return {
            'name': self.name,
            'subproject': self.subproject,
            'machine': self.machine,
        }

    def __setstate__(self, state: T.Dict[str, T.Any]) -> None:
        """De-serialize the state of a pickle.

        This is very clever. __init__ is not a constructor, it's an
        initializer, therefore it's safe to call more than once. We create a
        state in the custom __getstate__ method, which is valid to pass
        splatted to the initializer.
        """
        # Mypy doesn't like this, because it's so clever.
        self.__init__(**state)  # type: ignore

    def __hash__(self) -> int:
        return self._hash

    def _to_tuple(self) -> T.Tuple[str, str, str, MachineChoice, str]:
        return (self.subproject, self.machine, self.name)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, OptionKey):
            return self._to_tuple() == other._to_tuple()
        return NotImplemented

    def __lt__(self, other: object) -> bool:
        if isinstance(other, OptionKey):
            return self._to_tuple() < other._to_tuple()
        return NotImplemented

    def __str__(self) -> str:
        out = self.name
        if self.machine is MachineChoice.BUILD:
            out = f'build.{out}'
        if self.subproject:
            out = f'{self.subproject}:{out}'
        return out

    def __repr__(self) -> str:
        return f'OptionKey({self.name!r}, {self.subproject!r}, {self.machine!r})'

    @classmethod
    def from_string(cls, raw: str) -> 'OptionKey':
        """Parse the raw command line format into a three part tuple.

        This takes strings like `mysubproject:build.myoption` and Creates an
        OptionKey out of them.
        """
        try:
            subproject, raw2 = raw.split(':')
        except ValueError:
            subproject, raw2 = '', raw

        for_machine = MachineChoice.HOST
        try:
            prefix, raw3 = raw2.split('.')
            if prefix == 'build':
                for_machine = MachineChoice.BUILD
            else:
                raw3 = raw2
        except ValueError:
            raw3 = raw2

        opt = raw3
        assert ':' not in opt
        assert opt.count('.') < 2

        return cls(opt, subproject, for_machine)

    def evolve(self, name: T.Optional[str] = None, subproject: T.Optional[str] = None,
               machine: T.Optional[MachineChoice] = None) -> 'OptionKey':
        """Create a new copy of this key, but with altered members.

        For example:
        >>> a = OptionKey('foo', '', MachineChoice.Host)
        >>> b = OptionKey('foo', 'bar', MachineChoice.Host)
        >>> b == a.evolve(subproject='bar')
        True
        """
        # We have to be a little clever with lang here, because lang is valid
        # as None, for non-compiler options
        return OptionKey(
            name if name is not None else self.name,
            subproject if subproject is not None else self.subproject,
            machine if machine is not None else self.machine,
        )

    def as_root(self) -> 'OptionKey':
        """Convenience method for key.evolve(subproject='')."""
        return self.evolve(subproject='')

    def as_build(self) -> 'OptionKey':
        """Convenience method for key.evolve(machine=MachineChoice.BUILD)."""
        return self.evolve(machine=MachineChoice.BUILD)

    def as_host(self) -> 'OptionKey':
        """Convenience method for key.evolve(machine=MachineChoice.HOST)."""
        return self.evolve(machine=MachineChoice.HOST)

    def is_project_hack_for_optionsview(self) -> bool:
        """This method will be removed once we can delete OptionsView."""
        import sys
        sys.exit('FATAL internal error. This should not make it into an actual release. File a bug.')

    def has_module_prefix(self) -> bool:
        return '.' in self.name

    def get_module_prefix(self) -> T.Optional[str]:
        if self.has_module_prefix():
            return self.name.split('.', 1)[0]
        return None

    def without_module_prefix(self) -> 'OptionKey':
        if self.has_module_prefix():
            newname = self.name.split('.', 1)[1]
            return self.evolve(newname)
        return self


class UserOption(T.Generic[_T], HoldableObject):
    def __init__(self, name: str, description: str, choices: T.Optional[T.Union[str, T.List[_T]]],
                 yielding: bool,
                 deprecated: T.Union[bool, str, T.Dict[str, str], T.List[str]] = False):
        super().__init__()
        self.name = name
        self.choices = choices
        self.description = description
        if not isinstance(yielding, bool):
            raise MesonException('Value of "yielding" must be a boolean.')
        self.yielding = yielding
        self.deprecated = deprecated
        self.readonly = False

    def listify(self, value: T.Any) -> T.List[T.Any]:
        return [value]

    def printable_value(self) -> T.Union[str, int, bool, T.List[T.Union[str, int, bool]]]:
        assert isinstance(self.value, (str, int, bool, list))
        return self.value

    # Check that the input is a valid value and return the
    # "cleaned" or "native" version. For example the Boolean
    # option could take the string "true" and return True.
    def validate_value(self, value: T.Any) -> _T:
        raise RuntimeError('Derived option class did not override validate_value.')

    def set_value(self, newvalue: T.Any) -> bool:
        oldvalue = getattr(self, 'value', None)
        self.value = self.validate_value(newvalue)
        return self.value != oldvalue

_U = T.TypeVar('_U', bound=UserOption[_T])


class UserStringOption(UserOption[str]):
    def __init__(self, name: str, description: str, value: T.Any, yielding: bool = DEFAULT_YIELDING,
                 deprecated: T.Union[bool, str, T.Dict[str, str], T.List[str]] = False):
        super().__init__(name, description, None, yielding, deprecated)
        self.set_value(value)

    def validate_value(self, value: T.Any) -> str:
        if not isinstance(value, str):
            raise MesonException(f'The value of option "{self.name}" is "{value}", which is not a string.')
        return value

class UserBooleanOption(UserOption[bool]):
    def __init__(self, name: str, description: str, value: bool, yielding: bool = DEFAULT_YIELDING,
                 deprecated: T.Union[bool, str, T.Dict[str, str], T.List[str]] = False):
        super().__init__(name, description, [True, False], yielding, deprecated)
        self.set_value(value)

    def __bool__(self) -> bool:
        return self.value

    def validate_value(self, value: T.Any) -> bool:
        if isinstance(value, bool):
            return value
        if not isinstance(value, str):
            raise MesonException(f'Option "{self.name}" value {value} cannot be converted to a boolean')
        if value.lower() == 'true':
            return True
        if value.lower() == 'false':
            return False
        raise MesonException(f'Option "{self.name}" value {value} is not boolean (true or false).')

class UserIntegerOption(UserOption[int]):
    def __init__(self, name: str, description: str, value: T.Any, yielding: bool = DEFAULT_YIELDING,
                 deprecated: T.Union[bool, str, T.Dict[str, str], T.List[str]] = False):
        min_value, max_value, default_value = value
        self.min_value = min_value
        self.max_value = max_value
        c: T.List[str] = []
        if min_value is not None:
            c.append('>=' + str(min_value))
        if max_value is not None:
            c.append('<=' + str(max_value))
        choices = ', '.join(c)
        super().__init__(name, description, choices, yielding, deprecated)
        self.set_value(default_value)

    def validate_value(self, value: T.Any) -> int:
        if isinstance(value, str):
            value = self.toint(value)
        if not isinstance(value, int):
            raise MesonException(f'Value {value!r} for option "{self.name}" is not an integer.')
        if self.min_value is not None and value < self.min_value:
            raise MesonException(f'Value {value} for option "{self.name}" is less than minimum value {self.min_value}.')
        if self.max_value is not None and value > self.max_value:
            raise MesonException(f'Value {value} for option "{self.name}" is more than maximum value {self.max_value}.')
        return value

    def toint(self, valuestring: str) -> int:
        try:
            return int(valuestring)
        except ValueError:
            raise MesonException(f'Value string "{valuestring}" for option "{self.name}" is not convertible to an integer.')

class OctalInt(int):
    # NinjaBackend.get_user_option_args uses str() to converts it to a command line option
    # UserUmaskOption.toint() uses int(str, 8) to convert it to an integer
    # So we need to use oct instead of dec here if we do not want values to be misinterpreted.
    def __str__(self) -> str:
        return oct(int(self))

class UserUmaskOption(UserIntegerOption, UserOption[T.Union[str, OctalInt]]):
    def __init__(self, name: str, description: str, value: T.Any, yielding: bool = DEFAULT_YIELDING,
                 deprecated: T.Union[bool, str, T.Dict[str, str], T.List[str]] = False):
        super().__init__(name, description, (0, 0o777, value), yielding, deprecated)
        self.choices = ['preserve', '0000-0777']

    def printable_value(self) -> str:
        if self.value == 'preserve':
            return self.value
        return format(self.value, '04o')

    def validate_value(self, value: T.Any) -> T.Union[str, OctalInt]:
        if value == 'preserve':
            return 'preserve'
        return OctalInt(super().validate_value(value))

    def toint(self, valuestring: T.Union[str, OctalInt]) -> int:
        try:
            return int(valuestring, 8)
        except ValueError as e:
            raise MesonException(f'Invalid mode for option "{self.name}" {e}')

class UserComboOption(UserOption[str]):
    def __init__(self, name: str, description: str, choices: T.List[str], value: T.Any,
                 yielding: bool = DEFAULT_YIELDING,
                 deprecated: T.Union[bool, str, T.Dict[str, str], T.List[str]] = False):
        super().__init__(name, description, choices, yielding, deprecated)
        if not isinstance(self.choices, list):
            raise MesonException(f'Combo choices for option "{self.name}" must be an array.')
        for i in self.choices:
            if not isinstance(i, str):
                raise MesonException(f'Combo choice elements for option "{self.name}" must be strings.')
        self.set_value(value)

    def validate_value(self, value: T.Any) -> str:
        if value not in self.choices:
            if isinstance(value, bool):
                _type = 'boolean'
            elif isinstance(value, (int, float)):
                _type = 'number'
            else:
                _type = 'string'
            optionsstring = ', '.join([f'"{item}"' for item in self.choices])
            raise MesonException('Value "{}" (of type "{}") for option "{}" is not one of the choices.'
                                 ' Possible choices are (as string): {}.'.format(
                                     value, _type, self.name, optionsstring))
        return value

class UserArrayOption(UserOption[T.List[str]]):
    def __init__(self, name: str, description: str, value: T.Union[str, T.List[str]],
                 split_args: bool = False,
                 allow_dups: bool = False, yielding: bool = DEFAULT_YIELDING,
                 choices: T.Optional[T.List[str]] = None,
                 deprecated: T.Union[bool, str, T.Dict[str, str], T.List[str]] = False):
        super().__init__(name, description, choices if choices is not None else [], yielding, deprecated)
        self.split_args = split_args
        self.allow_dups = allow_dups
        self.set_value(value)

    def listify(self, value: T.Any) -> T.List[T.Any]:
        try:
            return listify_array_value(value, self.split_args)
        except MesonException as e:
            raise MesonException(f'error in option "{self.name}": {e!s}')

    def validate_value(self, value: T.Union[str, T.List[str]]) -> T.List[str]:
        newvalue = self.listify(value)

        if not self.allow_dups and len(set(newvalue)) != len(newvalue):
            msg = 'Duplicated values in array option is deprecated. ' \
                  'This will become a hard error in meson 2.0.'
            mlog.deprecation(msg)
        for i in newvalue:
            if not isinstance(i, str):
                raise MesonException(f'String array element "{newvalue!s}" for option "{self.name}" is not a string.')
        if self.choices:
            bad = [x for x in newvalue if x not in self.choices]
            if bad:
                raise MesonException('Value{} "{}" for option "{}" {} not in allowed choices: "{}"'.format(
                    '' if len(bad) == 1 else 's',
                    ', '.join(bad),
                    self.name,
                    'is' if len(bad) == 1 else 'are',
                    ', '.join(self.choices))
                )
        return newvalue

    def extend_value(self, value: T.Union[str, T.List[str]]) -> None:
        """Extend the value with an additional value."""
        new = self.validate_value(value)
        self.set_value(self.value + new)


class UserFeatureOption(UserComboOption):
    static_choices = ['enabled', 'disabled', 'auto']

    def __init__(self, name: str, description: str, value: T.Any, yielding: bool = DEFAULT_YIELDING,
                 deprecated: T.Union[bool, str, T.Dict[str, str], T.List[str]] = False):
        super().__init__(name, description, self.static_choices, value, yielding, deprecated)
        self.name: T.Optional[str] = None  # TODO: Refactor options to all store their name

    def is_enabled(self) -> bool:
        return self.value == 'enabled'

    def is_disabled(self) -> bool:
        return self.value == 'disabled'

    def is_auto(self) -> bool:
        return self.value == 'auto'

class UserStdOption(UserComboOption):
    '''
    UserOption specific to c_std and cpp_std options. User can set a list of
    STDs in preference order and it selects the first one supported by current
    compiler.

    For historical reasons, some compilers (msvc) allowed setting a GNU std and
    silently fell back to C std. This is now deprecated. Projects that support
    both GNU and MSVC compilers should set e.g. c_std=gnu11,c11.

    This is not using self.deprecated mechanism we already have for project
    options because we want to print a warning if ALL values are deprecated, not
    if SOME values are deprecated.
    '''
    def __init__(self, lang: str, all_stds: T.List[str]) -> None:
        self.lang = lang.lower()
        self.all_stds = ['none'] + all_stds
        # Map a deprecated std to its replacement. e.g. gnu11 -> c11.
        self.deprecated_stds: T.Dict[str, str] = {}
        opt_name = 'cpp_std' if lang == 'c++' else f'{lang}_std'
        super().__init__(opt_name, f'{lang} language standard to use', ['none'], 'none')

    def set_versions(self, versions: T.List[str], gnu: bool = False, gnu_deprecated: bool = False) -> None:
        assert all(std in self.all_stds for std in versions)
        self.choices += versions
        if gnu:
            gnu_stds_map = {f'gnu{std[1:]}': std for std in versions}
            if gnu_deprecated:
                self.deprecated_stds.update(gnu_stds_map)
            else:
                self.choices += gnu_stds_map.keys()

    def validate_value(self, value: T.Union[str, T.List[str]]) -> str:
        try:
            candidates = listify_array_value(value)
        except MesonException as e:
            raise MesonException(f'error in option "{self.name}": {e!s}')
        unknown = ','.join(std for std in candidates if std not in self.all_stds)
        if unknown:
            raise MesonException(f'Unknown option "{self.name}" value {unknown}. Possible values are {self.all_stds}.')
        # Check first if any of the candidates are not deprecated
        for std in candidates:
            if std in self.choices:
                return std
        # Fallback to a deprecated std if any
        for std in candidates:
            newstd = self.deprecated_stds.get(std)
            if newstd is not None:
                mlog.deprecation(
                    f'None of the values {candidates} are supported by the {self.lang} compiler.\n' +
                    f'However, the deprecated {std} std currently falls back to {newstd}.\n' +
                    'This will be an error in meson 2.0.\n' +
                    'If the project supports both GNU and MSVC compilers, a value such as\n' +
                    '"c_std=gnu11,c11" specifies that GNU is preferred but it can safely fallback to plain c11.')
                return newstd
        raise MesonException(f'None of values {candidates} are supported by the {self.lang.upper()} compiler. ' +
                             f'Possible values for option "{self.name}" are {self.choices}')


class BuiltinOption(T.Generic[_T, _U]):

    """Class for a builtin option type.

    There are some cases that are not fully supported yet.
    """

    def __init__(self, opt_type: T.Type[_U], description: str, default: T.Any, yielding: bool = True, *,
                 choices: T.Any = None, readonly: bool = False):
        self.opt_type = opt_type
        self.description = description
        self.default = default
        self.choices = choices
        self.yielding = yielding
        self.readonly = readonly

    def init_option(self, name: 'OptionKey', value: T.Optional[T.Any], prefix: str) -> _U:
        """Create an instance of opt_type and return it."""
        if value is None:
            value = self.prefixed_default(name, prefix)
        keywords = {'yielding': self.yielding, 'value': value}
        if self.choices:
            keywords['choices'] = self.choices
        o = self.opt_type(name.name, self.description, **keywords)
        o.readonly = self.readonly
        return o

    def _argparse_action(self) -> T.Optional[str]:
        # If the type is a boolean, the presence of the argument in --foo form
        # is to enable it. Disabling happens by using -Dfoo=false, which is
        # parsed under `args.projectoptions` and does not hit this codepath.
        if isinstance(self.default, bool):
            return 'store_true'
        return None

    def _argparse_choices(self) -> T.Any:
        if self.opt_type is UserBooleanOption:
            return [True, False]
        elif self.opt_type is UserFeatureOption:
            return UserFeatureOption.static_choices
        return self.choices

    @staticmethod
    def argparse_name_to_arg(name: str) -> str:
        if name == 'warning_level':
            return '--warnlevel'
        else:
            return '--' + name.replace('_', '-')

    def prefixed_default(self, name: 'OptionKey', prefix: str = '') -> T.Any:
        if self.opt_type in [UserComboOption, UserIntegerOption]:
            return self.default
        try:
            return BUILTIN_DIR_NOPREFIX_OPTIONS[name][prefix]
        except KeyError:
            pass
        return self.default

    def add_to_argparse(self, name: str, parser: argparse.ArgumentParser, help_suffix: str) -> None:
        kwargs: ArgparseKWs = {}

        c = self._argparse_choices()
        b = self._argparse_action()
        h = self.description
        if not b:
            h = '{} (default: {}).'.format(h.rstrip('.'), self.prefixed_default(name))
        else:
            kwargs['action'] = b
        if c and not b:
            kwargs['choices'] = c
        kwargs['default'] = argparse.SUPPRESS
        kwargs['dest'] = name

        cmdline_name = self.argparse_name_to_arg(name)
        parser.add_argument(cmdline_name, help=h + help_suffix, **kwargs)


# Update `docs/markdown/Builtin-options.md` after changing the options below
# Also update mesonlib._BUILTIN_NAMES. See the comment there for why this is required.
# Please also update completion scripts in $MESONSRC/data/shell-completions/
BUILTIN_DIR_OPTIONS: T.Dict['OptionKey', 'BuiltinOption'] = OrderedDict([
    (OptionKey('prefix'),          BuiltinOption(UserStringOption, 'Installation prefix', default_prefix())),
    (OptionKey('bindir'),          BuiltinOption(UserStringOption, 'Executable directory', 'bin')),
    (OptionKey('datadir'),         BuiltinOption(UserStringOption, 'Data file directory', default_datadir())),
    (OptionKey('includedir'),      BuiltinOption(UserStringOption, 'Header file directory', default_includedir())),
    (OptionKey('infodir'),         BuiltinOption(UserStringOption, 'Info page directory', default_infodir())),
    (OptionKey('libdir'),          BuiltinOption(UserStringOption, 'Library directory', default_libdir())),
    (OptionKey('licensedir'),      BuiltinOption(UserStringOption, 'Licenses directory', '')),
    (OptionKey('libexecdir'),      BuiltinOption(UserStringOption, 'Library executable directory', default_libexecdir())),
    (OptionKey('localedir'),       BuiltinOption(UserStringOption, 'Locale data directory', default_localedir())),
    (OptionKey('localstatedir'),   BuiltinOption(UserStringOption, 'Localstate data directory', 'var')),
    (OptionKey('mandir'),          BuiltinOption(UserStringOption, 'Manual page directory', default_mandir())),
    (OptionKey('sbindir'),         BuiltinOption(UserStringOption, 'System executable directory', default_sbindir())),
    (OptionKey('sharedstatedir'),  BuiltinOption(UserStringOption, 'Architecture-independent data directory', 'com')),
    (OptionKey('sysconfdir'),      BuiltinOption(UserStringOption, 'Sysconf data directory', default_sysconfdir())),
])

BUILTIN_CORE_OPTIONS: T.Dict['OptionKey', 'BuiltinOption'] = OrderedDict([
    (OptionKey('auto_features'),   BuiltinOption(UserFeatureOption, "Override value of all 'auto' features", 'auto')),
    (OptionKey('backend'),         BuiltinOption(UserComboOption, 'Backend to use', 'ninja', choices=backendlist,
                                                 readonly=True)),
    (OptionKey('genvslite'),
     BuiltinOption(
         UserComboOption,
         'Setup multiple buildtype-suffixed ninja-backend build directories, '
         'and a [builddir]_vs containing a Visual Studio meta-backend with multiple configurations that calls into them',
         'vs2022',
         choices=genvslitelist)
     ),
    (OptionKey('buildtype'),       BuiltinOption(UserComboOption, 'Build type to use', 'debug',
                                                 choices=buildtypelist)),
    (OptionKey('debug'),           BuiltinOption(UserBooleanOption, 'Enable debug symbols and other information', True)),
    (OptionKey('default_library'), BuiltinOption(UserComboOption, 'Default library type', 'shared', choices=['shared', 'static', 'both'],
                                                 yielding=False)),
    (OptionKey('default_both_libraries'), BuiltinOption(UserComboOption, 'Default library type for both_libraries', 'shared', choices=['shared', 'static', 'auto'])),
    (OptionKey('errorlogs'),       BuiltinOption(UserBooleanOption, "Whether to print the logs from failing tests", True)),
    (OptionKey('install_umask'),   BuiltinOption(UserUmaskOption, 'Default umask to apply on permissions of installed files', '022')),
    (OptionKey('layout'),          BuiltinOption(UserComboOption, 'Build directory layout', 'mirror', choices=['mirror', 'flat'])),
    (OptionKey('optimization'),    BuiltinOption(UserComboOption, 'Optimization level', '0', choices=['plain', '0', 'g', '1', '2', '3', 's'])),
    (OptionKey('prefer_static'),   BuiltinOption(UserBooleanOption, 'Whether to try static linking before shared linking', False)),
    (OptionKey('stdsplit'),        BuiltinOption(UserBooleanOption, 'Split stdout and stderr in test logs', True)),
    (OptionKey('strip'),           BuiltinOption(UserBooleanOption, 'Strip targets on install', False)),
    (OptionKey('unity'),           BuiltinOption(UserComboOption, 'Unity build', 'off', choices=['on', 'off', 'subprojects'])),
    (OptionKey('unity_size'),      BuiltinOption(UserIntegerOption, 'Unity block size', (2, None, 4))),
    (OptionKey('warning_level'),   BuiltinOption(UserComboOption, 'Compiler warning level to use', '1', choices=['0', '1', '2', '3', 'everything'], yielding=False)),
    (OptionKey('werror'),          BuiltinOption(UserBooleanOption, 'Treat warnings as errors', False, yielding=False)),
    (OptionKey('wrap_mode'),       BuiltinOption(UserComboOption, 'Wrap mode', 'default', choices=['default', 'nofallback', 'nodownload', 'forcefallback', 'nopromote'])),
    (OptionKey('force_fallback_for'), BuiltinOption(UserArrayOption, 'Force fallback for those subprojects', [])),
    (OptionKey('vsenv'),           BuiltinOption(UserBooleanOption, 'Activate Visual Studio environment', False, readonly=True)),

    # Pkgconfig module
    (OptionKey('pkgconfig.relocatable'),
     BuiltinOption(UserBooleanOption, 'Generate pkgconfig files as relocatable', False)),

    # Python module
    (OptionKey('python.bytecompile'),
     BuiltinOption(UserIntegerOption, 'Whether to compile bytecode', (-1, 2, 0))),
    (OptionKey('python.install_env'),
     BuiltinOption(UserComboOption, 'Which python environment to install to', 'prefix', choices=['auto', 'prefix', 'system', 'venv'])),
    (OptionKey('python.platlibdir'),
     BuiltinOption(UserStringOption, 'Directory for site-specific, platform-specific files.', '')),
    (OptionKey('python.purelibdir'),
     BuiltinOption(UserStringOption, 'Directory for site-specific, non-platform-specific files.', '')),
    (OptionKey('python.allow_limited_api'),
     BuiltinOption(UserBooleanOption, 'Whether to allow use of the Python Limited API', True)),
])

BUILTIN_OPTIONS = OrderedDict(chain(BUILTIN_DIR_OPTIONS.items(), BUILTIN_CORE_OPTIONS.items()))

BUILTIN_OPTIONS_PER_MACHINE: T.Dict['OptionKey', 'BuiltinOption'] = OrderedDict([
    (OptionKey('pkg_config_path'), BuiltinOption(UserArrayOption, 'List of additional paths for pkg-config to search', [])),
    (OptionKey('cmake_prefix_path'), BuiltinOption(UserArrayOption, 'List of additional prefixes for cmake to search', [])),
])

# Special prefix-dependent defaults for installation directories that reside in
# a path outside of the prefix in FHS and common usage.
BUILTIN_DIR_NOPREFIX_OPTIONS: T.Dict[OptionKey, T.Dict[str, str]] = {
    OptionKey('sysconfdir'):     {'/usr': '/etc'},
    OptionKey('localstatedir'):  {'/usr': '/var',     '/usr/local': '/var/local'},
    OptionKey('sharedstatedir'): {'/usr': '/var/lib', '/usr/local': '/var/local/lib'},
    OptionKey('python.platlibdir'): {},
    OptionKey('python.purelibdir'): {},
}

class OptionStore:
    def __init__(self) -> None:
        self.d: T.Dict['OptionKey', 'UserOption[T.Any]'] = {}
        self.project_options: T.Set[OptionKey] = set()
        self.module_options: T.Set[OptionKey] = set()
        from .compilers import all_languages
        self.all_languages = set(all_languages)

    def __len__(self) -> int:
        return len(self.d)

    def ensure_key(self, key: T.Union[OptionKey, str]) -> OptionKey:
        if isinstance(key, str):
            return OptionKey(key)
        return key

    def get_value_object(self, key: T.Union[OptionKey, str]) -> 'UserOption[T.Any]':
        return self.d[self.ensure_key(key)]

    def get_value(self, key: T.Union[OptionKey, str]) -> 'T.Any':
        return self.get_value_object(key).value

    def add_system_option(self, key: T.Union[OptionKey, str], valobj: 'UserOption[T.Any]') -> None:
        key = self.ensure_key(key)
        if '.' in key.name:
            raise MesonException(f'Internal error: non-module option has a period in its name {key.name}.')
        self.add_system_option_internal(key, valobj)

    def add_system_option_internal(self, key: T.Union[OptionKey, str], valobj: 'UserOption[T.Any]') -> None:
        key = self.ensure_key(key)
        assert isinstance(valobj, UserOption)
        self.d[key] = valobj

    def add_compiler_option(self, language: str, key: T.Union[OptionKey, str], valobj: 'UserOption[T.Any]') -> None:
        key = self.ensure_key(key)
        if not key.name.startswith(language + '_'):
            raise MesonException(f'Internal error: all compiler option names must start with language prefix. ({key.name} vs {language}_)')
        self.add_system_option(key, valobj)

    def add_project_option(self, key: T.Union[OptionKey, str], valobj: 'UserOption[T.Any]') -> None:
        key = self.ensure_key(key)
        self.d[key] = valobj
        self.project_options.add(key)

    def add_module_option(self, modulename: str, key: T.Union[OptionKey, str], valobj: 'UserOption[T.Any]') -> None:
        key = self.ensure_key(key)
        if key.name.startswith('build.'):
            raise MesonException('FATAL internal error: somebody goofed option handling.')
        if not key.name.startswith(modulename + '.'):
            raise MesonException('Internal error: module option name {key.name} does not start with module prefix {modulename}.')
        self.add_system_option_internal(key, valobj)
        self.module_options.add(key)

    def set_value(self, key: T.Union[OptionKey, str], new_value: 'T.Any') -> bool:
        key = self.ensure_key(key)
        return self.d[key].set_value(new_value)

    # FIXME, this should be removed.or renamed to "change_type_of_existing_object" or something like that
    def set_value_object(self, key: T.Union[OptionKey, str], new_object: 'UserOption[T.Any]') -> None:
        key = self.ensure_key(key)
        self.d[key] = new_object

    def remove(self, key: OptionKey) -> None:
        del self.d[key]

    def __contains__(self, key: OptionKey) -> bool:
        key = self.ensure_key(key)
        return key in self.d

    def __repr__(self) -> str:
        return repr(self.d)

    def keys(self) -> T.KeysView[OptionKey]:
        return self.d.keys()

    def values(self) -> T.ValuesView[UserOption[T.Any]]:
        return self.d.values()

    def items(self) -> T.ItemsView['OptionKey', 'UserOption[T.Any]']:
        return self.d.items()

    # FIXME: this method must be deleted and users moved to use "add_xxx_option"s instead.
    def update(self, **kwargs: UserOption[T.Any]) -> None:
        self.d.update(**kwargs)

    def setdefault(self, k: OptionKey, o: UserOption[T.Any]) -> UserOption[T.Any]:
        return self.d.setdefault(k, o)

    def get(self, o: OptionKey, default: T.Optional[UserOption[T.Any]] = None) -> T.Optional[UserOption[T.Any]]:
        return self.d.get(o, default)

    def is_project_option(self, key: OptionKey) -> bool:
        """Convenience method to check if this is a project option."""
        return key in self.project_options

    def is_reserved_name(self, key: OptionKey) -> bool:
        if key.name in _BUILTIN_NAMES:
            return True
        if '_' not in key.name:
            return False
        prefix = key.name.split('_')[0]
        # Pylint seems to think that it is faster to build a set object
        # and all related work just to test whether a string has one of two
        # values. It is not, thank you very much.
        if prefix in ('b', 'backend'): # pylint: disable=R6201
            return True
        if prefix in self.all_languages:
            return True
        return False

    def is_builtin_option(self, key: OptionKey) -> bool:
        """Convenience method to check if this is a builtin option."""
        return key.name in _BUILTIN_NAMES or self.is_module_option(key)

    def is_base_option(self, key: OptionKey) -> bool:
        """Convenience method to check if this is a base option."""
        return key.name.startswith('b_')

    def is_backend_option(self, key: OptionKey) -> bool:
        """Convenience method to check if this is a backend option."""
        return key.name.startswith('backend_')

    def is_compiler_option(self, key: OptionKey) -> bool:
        """Convenience method to check if this is a compiler option."""

        # FIXME, duplicate of is_reserved_name above. Should maybe store a cache instead.
        if '_' not in key.name:
            return False
        prefix = key.name.split('_')[0]
        if prefix in self.all_languages:
            return True
        return False

    def is_module_option(self, key: OptionKey) -> bool:
        return key in self.module_options
