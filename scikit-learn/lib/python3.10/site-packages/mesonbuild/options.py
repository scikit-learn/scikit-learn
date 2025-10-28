# SPDX-License-Identifier: Apache-2.0
# Copyright 2013-2024 Contributors to the The Meson project
# Copyright © 2019-2025 Intel Corporation

from __future__ import annotations
from collections import OrderedDict
from itertools import chain
import argparse
import copy
import dataclasses
import itertools
import os
import pathlib

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
    MesonBugException,
    listify_array_value,
    MachineChoice,
)
from . import mlog

if T.TYPE_CHECKING:
    from typing_extensions import Literal, Final, TypeAlias, TypedDict

    from .interpreterbase import SubProject

    DeprecatedType: TypeAlias = T.Union[bool, str, T.Dict[str, str], T.List[str]]
    AnyOptionType: TypeAlias = T.Union[
        'UserBooleanOption', 'UserComboOption', 'UserFeatureOption',
        'UserIntegerOption', 'UserStdOption', 'UserStringArrayOption',
        'UserStringOption', 'UserUmaskOption']
    ElementaryOptionValues: TypeAlias = T.Union[str, int, bool, T.List[str]]
    MutableKeyedOptionDictType: TypeAlias = T.Dict['OptionKey', AnyOptionType]

    _OptionKeyTuple: TypeAlias = T.Tuple[T.Optional[str], MachineChoice, str]

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

_BAD_VALUE = 'Qwert Zuiopü'
_optionkey_cache: T.Dict[_OptionKeyTuple, OptionKey] = {}


class OptionKey:

    """Represents an option key in the various option dictionaries.

    This provides a flexible, powerful way to map option names from their
    external form (things like subproject:build.option) to something that
    internally easier to reason about and produce.
    """

    __slots__ = ('name', 'subproject', 'machine', '_hash')

    name: str
    subproject: T.Optional[str]  # None is global, empty string means top level project
    machine: MachineChoice
    _hash: int

    def __new__(cls,
                name: str = '',
                subproject: T.Optional[str] = None,
                machine: MachineChoice = MachineChoice.HOST) -> OptionKey:
        """The use of the __new__ method allows to add a transparent cache
        to the OptionKey object creation, without breaking its API.
        """
        if not name:
            return super().__new__(cls)  # for unpickling, do not cache now

        tuple_: _OptionKeyTuple = (subproject, machine, name)
        try:
            return _optionkey_cache[tuple_]
        except KeyError:
            instance = super().__new__(cls)
            instance._init(name, subproject, machine)
            _optionkey_cache[tuple_] = instance
            return instance

    def _init(self, name: str, subproject: T.Optional[str], machine: MachineChoice) -> None:
        # We don't use the __init__ method, because it would be called after __new__
        # while we need __new__ to initialise the object before populating the cache.

        if not isinstance(machine, MachineChoice):
            raise MesonException(f'Internal error, bad machine type: {machine}')
        if not isinstance(name, str):
            raise MesonBugException(f'Key name is not a string: {name}')
        assert ':' not in name

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
        # Here, the object is created using __new__()
        self._init(**state)
        _optionkey_cache[self._to_tuple()] = self

    def __hash__(self) -> int:
        return self._hash

    def _to_tuple(self) -> _OptionKeyTuple:
        return (self.subproject, self.machine, self.name)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, OptionKey):
            return self._to_tuple() == other._to_tuple()
        return NotImplemented

    def __ne__(self, other: object) -> bool:
        if isinstance(other, OptionKey):
            return self._to_tuple() != other._to_tuple()
        return NotImplemented

    def __lt__(self, other: object) -> bool:
        if isinstance(other, OptionKey):
            if self.subproject is None:
                return other.subproject is not None
            elif other.subproject is None:
                return False
            return self._to_tuple() < other._to_tuple()
        return NotImplemented

    def __le__(self, other: object) -> bool:
        if isinstance(other, OptionKey):
            if self.subproject is None and other.subproject is not None:
                return True
            elif self.subproject is not None and other.subproject is None:
                return False
            return self._to_tuple() <= other._to_tuple()
        return NotImplemented

    def __gt__(self, other: object) -> bool:
        if isinstance(other, OptionKey):
            if other.subproject is None:
                return self.subproject is not None
            elif self.subproject is None:
                return False
            return self._to_tuple() > other._to_tuple()
        return NotImplemented

    def __ge__(self, other: object) -> bool:
        if isinstance(other, OptionKey):
            if self.subproject is None and other.subproject is not None:
                return False
            elif self.subproject is not None and other.subproject is None:
                return True
            return self._to_tuple() >= other._to_tuple()
        return NotImplemented

    def __str__(self) -> str:
        out = self.name
        if self.machine is MachineChoice.BUILD:
            out = f'build.{out}'
        if self.subproject is not None:
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
        assert isinstance(raw, str)
        try:
            subproject, raw2 = raw.split(':')
        except ValueError:
            subproject, raw2 = None, raw

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

    def evolve(self,
               name: T.Optional[str] = None,
               subproject: T.Optional[str] = _BAD_VALUE,
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
        return OptionKey(name if name is not None else self.name,
                         subproject if subproject != _BAD_VALUE else self.subproject, # None is a valid value so it can'the default value in method declaration.
                         machine if machine is not None else self.machine)

    def as_root(self) -> OptionKey:
        """Convenience method for key.evolve(subproject='')."""
        return self.evolve(subproject='')

    def as_build(self) -> OptionKey:
        """Convenience method for key.evolve(machine=MachineChoice.BUILD)."""
        return self.evolve(machine=MachineChoice.BUILD)

    def as_host(self) -> OptionKey:
        """Convenience method for key.evolve(machine=MachineChoice.HOST)."""
        return self.evolve(machine=MachineChoice.HOST)

    def has_module_prefix(self) -> bool:
        return '.' in self.name

    def get_module_prefix(self) -> T.Optional[str]:
        if self.has_module_prefix():
            return self.name.split('.', 1)[0]
        return None

    def is_for_build(self) -> bool:
        return self.machine is MachineChoice.BUILD

if T.TYPE_CHECKING:
    OptionDict: TypeAlias = T.Dict[OptionKey, ElementaryOptionValues]

@dataclasses.dataclass
class UserOption(T.Generic[_T], HoldableObject):

    name: str
    description: str
    value_: dataclasses.InitVar[_T]
    yielding: bool = DEFAULT_YIELDING
    deprecated: DeprecatedType = False
    readonly: bool = dataclasses.field(default=False)
    parent: T.Optional[UserOption] = None

    def __post_init__(self, value_: _T) -> None:
        self.value = self.validate_value(value_)
        # Final isn't technically allowed in a __post_init__ method
        self.default: Final[_T] = self.value  # type: ignore[misc]

    def listify(self, value: ElementaryOptionValues) -> T.List[str]:
        if isinstance(value, list):
            return value
        if isinstance(value, bool):
            return ['true'] if value else ['false']
        if isinstance(value, int):
            return [str(value)]
        return [value]

    def printable_value(self) -> ElementaryOptionValues:
        assert isinstance(self.value, (str, int, bool, list))
        return self.value

    def printable_choices(self) -> T.Optional[T.List[str]]:
        return None

    # Check that the input is a valid value and return the
    # "cleaned" or "native" version. For example the Boolean
    # option could take the string "true" and return True.
    def validate_value(self, value: object) -> _T:
        raise RuntimeError('Derived option class did not override validate_value.')

    def set_value(self, newvalue: object) -> bool:
        oldvalue = self.value
        self.value = self.validate_value(newvalue)
        return self.value != oldvalue

@dataclasses.dataclass
class EnumeratedUserOption(UserOption[_T]):

    """A generic UserOption that has enumerated values."""

    choices: T.List[_T] = dataclasses.field(default_factory=list)

    def printable_choices(self) -> T.Optional[T.List[str]]:
        return [str(c) for c in self.choices]


class UserStringOption(UserOption[str]):

    def validate_value(self, value: object) -> str:
        if not isinstance(value, str):
            raise MesonException(f'The value of option "{self.name}" is "{value}", which is not a string.')
        return value

@dataclasses.dataclass
class UserBooleanOption(EnumeratedUserOption[bool]):

    choices: T.List[bool] = dataclasses.field(default_factory=lambda: [True, False])

    def __bool__(self) -> bool:
        return self.value

    def validate_value(self, value: object) -> bool:
        if isinstance(value, bool):
            return value
        if not isinstance(value, str):
            raise MesonException(f'Option "{self.name}" value {value} cannot be converted to a boolean')
        if value.lower() == 'true':
            return True
        if value.lower() == 'false':
            return False
        raise MesonException(f'Option "{self.name}" value {value} is not boolean (true or false).')


class _UserIntegerBase(UserOption[_T]):

    min_value: T.Optional[int]
    max_value: T.Optional[int]

    if T.TYPE_CHECKING:
        def toint(self, v: str) -> int: ...

    def __post_init__(self, value_: _T) -> None:
        super().__post_init__(value_)
        choices: T.List[str] = []
        if self.min_value is not None:
            choices.append(f'>= {self.min_value!s}')
        if self.max_value is not None:
            choices.append(f'<= {self.max_value!s}')
        self.__choices: str = ', '.join(choices)

    def printable_choices(self) -> T.Optional[T.List[str]]:
        return [self.__choices]

    def validate_value(self, value: object) -> _T:
        if isinstance(value, str):
            value = T.cast('_T', self.toint(value))
        if not isinstance(value, int):
            raise MesonException(f'Value {value!r} for option "{self.name}" is not an integer.')
        if self.min_value is not None and value < self.min_value:
            raise MesonException(f'Value {value} for option "{self.name}" is less than minimum value {self.min_value}.')
        if self.max_value is not None and value > self.max_value:
            raise MesonException(f'Value {value} for option "{self.name}" is more than maximum value {self.max_value}.')
        return T.cast('_T', value)


@dataclasses.dataclass
class UserIntegerOption(_UserIntegerBase[int]):

    min_value: T.Optional[int] = None
    max_value: T.Optional[int] = None

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


@dataclasses.dataclass
class UserUmaskOption(_UserIntegerBase[T.Union["Literal['preserve']", OctalInt]]):

    min_value: T.Optional[int] = dataclasses.field(default=0, init=False)
    max_value: T.Optional[int] = dataclasses.field(default=0o777, init=False)

    def printable_value(self) -> str:
        if isinstance(self.value, int):
            return format(self.value, '04o')
        return self.value

    def validate_value(self, value: object) -> T.Union[Literal['preserve'], OctalInt]:
        if value == 'preserve':
            return 'preserve'
        return OctalInt(super().validate_value(value))

    def toint(self, valuestring: str) -> int:
        try:
            return int(valuestring, 8)
        except ValueError as e:
            raise MesonException(f'Invalid mode for option "{self.name}" {e}')


@dataclasses.dataclass
class UserComboOption(EnumeratedUserOption[str]):

    def validate_value(self, value: object) -> str:
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

        assert isinstance(value, str), 'for mypy'
        return value

@dataclasses.dataclass
class UserArrayOption(UserOption[T.List[_T]]):

    value_: dataclasses.InitVar[T.Union[_T, T.List[_T]]]
    choices: T.Optional[T.List[_T]] = None
    split_args: bool = False
    allow_dups: bool = False

    def extend_value(self, value: T.Union[str, T.List[str]]) -> None:
        """Extend the value with an additional value."""
        new = self.validate_value(value)
        self.set_value(self.value + new)

    def printable_choices(self) -> T.Optional[T.List[str]]:
        if self.choices is None:
            return None
        return [str(c) for c in self.choices]


@dataclasses.dataclass
class UserStringArrayOption(UserArrayOption[str]):

    def listify(self, value: object) -> T.List[str]:
        try:
            return listify_array_value(value, self.split_args)
        except MesonException as e:
            raise MesonException(f'error in option "{self.name}": {e!s}')

    def validate_value(self, value: object) -> T.List[str]:
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


@dataclasses.dataclass
class UserFeatureOption(UserComboOption):

    choices: T.List[str] = dataclasses.field(
        # Ensure we get a copy with the lambda
        default_factory=lambda: ['enabled', 'disabled', 'auto'], init=False)

    def is_enabled(self) -> bool:
        return self.value == 'enabled'

    def is_disabled(self) -> bool:
        return self.value == 'disabled'

    def is_auto(self) -> bool:
        return self.value == 'auto'


_U = T.TypeVar('_U', bound=UserOption)


def choices_are_different(a: _U, b: _U) -> bool:
    """Are the choices between two options the same?

    :param a: A UserOption[T]
    :param b: A second UserOption[T]
    :return: True if the choices have changed, otherwise False
    """
    if isinstance(a, EnumeratedUserOption):
        # We expect `a` and `b` to be of the same type, but can't really annotate it that way.
        assert isinstance(b, EnumeratedUserOption), 'for mypy'
        return a.choices != b.choices
    elif isinstance(a, UserArrayOption):
        # We expect `a` and `b` to be of the same type, but can't really annotate it that way.
        assert isinstance(b, UserArrayOption), 'for mypy'
        return a.choices != b.choices
    elif isinstance(a, _UserIntegerBase):
        assert isinstance(b, _UserIntegerBase), 'for mypy'
        return a.max_value != b.max_value or a.min_value != b.min_value

    return False


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
        super().__init__(opt_name, f'{lang} language standard to use', 'none', choices=['none'])

    def set_versions(self, versions: T.List[str], gnu: bool = False, gnu_deprecated: bool = False) -> None:
        assert all(std in self.all_stds for std in versions)
        self.choices += versions
        if gnu:
            gnu_stds_map = {f'gnu{std[1:]}': std for std in versions}
            if gnu_deprecated:
                self.deprecated_stds.update(gnu_stds_map)
            else:
                self.choices += gnu_stds_map.keys()

    def validate_value(self, value: object) -> str:
        try:
            candidates = listify_array_value(value)
        except MesonException as e:
            raise MesonException(f'error in option "{self.name}": {e!s}')
        for std in candidates:
            if not isinstance(std, str):
                raise MesonException(f'String array element "{candidates!s}" for option "{self.name}" is not a string.')
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
                    '"c_std=gnu11,c11" specifies that GNU is preferred but it can safely fallback to plain c11.', once=True)
                return newstd
        raise MesonException(f'None of values {candidates} are supported by the {self.lang.upper()} compiler. ' +
                             f'Possible values for option "{self.name}" are {self.choices}')


def argparse_name_to_arg(name: str) -> str:
    if name == 'warning_level':
        return '--warnlevel'
    return '--' + name.replace('_', '-')


def argparse_prefixed_default(opt: AnyOptionType, name: OptionKey, prefix: str = '') -> ElementaryOptionValues:
    if isinstance(opt, (UserComboOption, UserIntegerOption, UserUmaskOption)):
        return T.cast('ElementaryOptionValues', opt.default)
    try:
        return BUILTIN_DIR_NOPREFIX_OPTIONS[name][prefix]
    except KeyError:
        return T.cast('ElementaryOptionValues', opt.default)


def option_to_argparse(option: AnyOptionType, name: OptionKey, parser: argparse.ArgumentParser, help_suffix: str) -> None:
    kwargs: ArgparseKWs = {}

    if isinstance(option, (EnumeratedUserOption, UserArrayOption)):
        c = option.choices
    else:
        c = None
    b = 'store_true' if isinstance(option.default, bool) else None
    h = option.description
    if not b:
        h = '{} (default: {}).'.format(h.rstrip('.'), argparse_prefixed_default(option, name))
    else:
        kwargs['action'] = b
    if c and not b:
        kwargs['choices'] = c
    kwargs['default'] = argparse.SUPPRESS
    kwargs['dest'] = str(name)

    cmdline_name = argparse_name_to_arg(str(name))
    parser.add_argument(cmdline_name, help=h + help_suffix, **kwargs)


# Update `docs/markdown/Builtin-options.md` after changing the options below
# Also update mesonlib._BUILTIN_NAMES. See the comment there for why this is required.
# Please also update completion scripts in $MESONSRC/data/shell-completions/
BUILTIN_DIR_OPTIONS: T.Mapping[OptionKey, AnyOptionType] = {
    OptionKey(o.name): o for o in [
        UserStringOption('prefix', 'Installation prefix', default_prefix()),
        UserStringOption('bindir', 'Executable directory', 'bin'),
        UserStringOption('datadir', 'Data file directory', default_datadir()),
        UserStringOption('includedir', 'Header file directory', default_includedir()),
        UserStringOption('infodir', 'Info page directory', default_infodir()),
        UserStringOption('libdir', 'Library directory', default_libdir()),
        UserStringOption('licensedir', 'Licenses directory', ''),
        UserStringOption('libexecdir', 'Library executable directory', default_libexecdir()),
        UserStringOption('localedir', 'Locale data directory', default_localedir()),
        UserStringOption('localstatedir', 'Localstate data directory', 'var'),
        UserStringOption('mandir', 'Manual page directory', default_mandir()),
        UserStringOption('sbindir', 'System executable directory', default_sbindir()),
        UserStringOption('sharedstatedir', 'Architecture-independent data directory', 'com'),
        UserStringOption('sysconfdir', 'Sysconf data directory', default_sysconfdir()),
    ]
}

BUILTIN_CORE_OPTIONS: T.Mapping[OptionKey, AnyOptionType] = {
    OptionKey(o.name): o for o in T.cast('T.List[AnyOptionType]', [
        UserFeatureOption('auto_features', "Override value of all 'auto' features", 'auto'),
        UserComboOption('backend', 'Backend to use', 'ninja', choices=backendlist, readonly=True),
        UserComboOption(
            'genvslite',
            'Setup multiple buildtype-suffixed ninja-backend build directories, '
            'and a [builddir]_vs containing a Visual Studio meta-backend with multiple configurations that calls into them',
            'vs2022',
            choices=genvslitelist
        ),
        UserComboOption('buildtype', 'Build type to use', 'debug', choices=buildtypelist),
        UserBooleanOption('debug', 'Enable debug symbols and other information', True),
        UserComboOption('default_library', 'Default library type', 'shared', choices=['shared', 'static', 'both'],
                        yielding=False),
        UserComboOption('default_both_libraries', 'Default library type for both_libraries', 'shared',
                        choices=['shared', 'static', 'auto']),
        UserBooleanOption('errorlogs', "Whether to print the logs from failing tests", True),
        UserUmaskOption('install_umask', 'Default umask to apply on permissions of installed files', OctalInt(0o022)),
        UserComboOption('layout', 'Build directory layout', 'mirror', choices=['mirror', 'flat']),
        UserComboOption('optimization', 'Optimization level', '0', choices=['plain', '0', 'g', '1', '2', '3', 's']),
        UserBooleanOption('prefer_static', 'Whether to try static linking before shared linking', False),
        UserBooleanOption('stdsplit', 'Split stdout and stderr in test logs', True),
        UserBooleanOption('strip', 'Strip targets on install', False),
        UserComboOption('unity', 'Unity build', 'off', choices=['on', 'off', 'subprojects']),
        UserIntegerOption('unity_size', 'Unity block size', 4, min_value=2),
        UserComboOption('warning_level', 'Compiler warning level to use', '1', choices=['0', '1', '2', '3', 'everything'],
                        yielding=False),
        UserBooleanOption('werror', 'Treat warnings as errors', False, yielding=False),
        UserComboOption('wrap_mode', 'Wrap mode', 'default', choices=['default', 'nofallback', 'nodownload', 'forcefallback', 'nopromote']),
        UserStringArrayOption('force_fallback_for', 'Force fallback for those subprojects', []),
        UserBooleanOption('vsenv', 'Activate Visual Studio environment', False, readonly=True),

        # Pkgconfig module
        UserBooleanOption('pkgconfig.relocatable', 'Generate pkgconfig files as relocatable', False),

        # Python module
        UserIntegerOption('python.bytecompile', 'Whether to compile bytecode', 0, min_value=-1, max_value=2),
        UserComboOption('python.install_env', 'Which python environment to install to', 'prefix',
                        choices=['auto', 'prefix', 'system', 'venv']),
        UserStringOption('python.platlibdir', 'Directory for site-specific, platform-specific files.', ''),
        UserStringOption('python.purelibdir', 'Directory for site-specific, non-platform-specific files.', ''),
        UserBooleanOption('python.allow_limited_api', 'Whether to allow use of the Python Limited API', True),
    ])
}

BUILTIN_OPTIONS = OrderedDict(chain(BUILTIN_DIR_OPTIONS.items(), BUILTIN_CORE_OPTIONS.items()))

BUILTIN_OPTIONS_PER_MACHINE: T.Mapping[OptionKey, AnyOptionType] = {
    OptionKey(o.name): o for o in [
        UserStringArrayOption('pkg_config_path', 'List of additional paths for pkg-config to search', []),
        UserStringArrayOption('cmake_prefix_path', 'List of additional prefixes for cmake to search', []),
    ]
}

# Special prefix-dependent defaults for installation directories that reside in
# a path outside of the prefix in FHS and common usage.
BUILTIN_DIR_NOPREFIX_OPTIONS: T.Dict[OptionKey, T.Dict[str, str]] = {
    OptionKey('sysconfdir'):     {'/usr': '/etc'},
    OptionKey('localstatedir'):  {'/usr': '/var',     '/usr/local': '/var/local'},
    OptionKey('sharedstatedir'): {'/usr': '/var/lib', '/usr/local': '/var/local/lib'},
    OptionKey('python.platlibdir'): {},
    OptionKey('python.purelibdir'): {},
}

MSCRT_VALS = ['none', 'md', 'mdd', 'mt', 'mtd']

COMPILER_BASE_OPTIONS: T.Mapping[OptionKey, AnyOptionType] = {
    OptionKey(o.name): o for o in T.cast('T.List[AnyOptionType]', [
        UserBooleanOption('b_pch', 'Use precompiled headers', True),
        UserBooleanOption('b_lto', 'Use link time optimization', False),
        UserIntegerOption('b_lto_threads', 'Use multiple threads for Link Time Optimization', 0),
        UserComboOption('b_lto_mode', 'Select between different LTO modes.', 'default', choices=['default', 'thin']),
        UserBooleanOption('b_thinlto_cache', 'Use LLVM ThinLTO caching for faster incremental builds', False),
        UserStringOption('b_thinlto_cache_dir', 'Directory to store ThinLTO cache objects', ''),
        UserStringArrayOption('b_sanitize', 'Code sanitizer to use', []),
        UserBooleanOption('b_lundef', 'Use -Wl,--no-undefined when linking', True),
        UserBooleanOption('b_asneeded', 'Use -Wl,--as-needed when linking', True),
        UserComboOption(
            'b_pgo', 'Use profile guided optimization', 'off', choices=['off', 'generate', 'use']),
        UserBooleanOption('b_coverage', 'Enable coverage tracking.', False),
        UserComboOption(
            'b_colorout', 'Use colored output', 'always', choices=['auto', 'always', 'never']),
        UserComboOption(
            'b_ndebug', 'Disable asserts', 'false', choices=['true', 'false', 'if-release']),
        UserBooleanOption('b_staticpic', 'Build static libraries as position independent', True),
        UserBooleanOption('b_pie', 'Build executables as position independent', False),
        UserBooleanOption('b_bitcode', 'Generate and embed bitcode (only macOS/iOS/tvOS)', False),
        UserComboOption(
            'b_vscrt', 'VS run-time library type to use.', 'from_buildtype',
            choices=MSCRT_VALS + ['from_buildtype', 'static_from_buildtype']),
    ])
}

class OptionStore:
    DEFAULT_DEPENDENTS = {'plain': ('plain', False),
                          'debug': ('0', True),
                          'debugoptimized': ('2', True),
                          'release': ('3', False),
                          'minsize': ('s', True),
                          }

    def __init__(self, is_cross: bool) -> None:
        self.options: T.Dict['OptionKey', 'AnyOptionType'] = {}
        self.subprojects: T.Set[str] = set()
        self.project_options: T.Set[OptionKey] = set()
        self.module_options: T.Set[OptionKey] = set()
        from .compilers import all_languages
        self.all_languages = set(all_languages)
        self.augments: OptionDict = {}
        self.is_cross = is_cross

        # Pending options are configuration dependent options that could be
        # initialized later, such as compiler options
        self.pending_options: OptionDict = {}
        # Subproject options from toplevel project()
        self.pending_subproject_options: OptionDict = {}

    def ensure_and_validate_key(self, key: T.Union[OptionKey, str]) -> OptionKey:
        if isinstance(key, str):
            return OptionKey(key)
        # FIXME. When not cross building all "build" options need to fall back
        # to "host" options due to how the old code worked.
        #
        # This is NOT how it should be.
        #
        # This needs to be changed to that trying to add or access "build" keys
        # is a hard error and fix issues that arise.
        #
        # I did not do this yet, because it would make this MR even
        # more massive than it already is. Later then.
        if not self.is_cross and key.machine == MachineChoice.BUILD:
            key = key.as_host()
        return key

    def get_pending_value(self, key: T.Union[OptionKey, str], default: T.Optional[ElementaryOptionValues] = None) -> ElementaryOptionValues:
        key = self.ensure_and_validate_key(key)
        if key in self.options:
            return self.options[key].value
        return self.pending_options.get(key, default)

    def get_value(self, key: T.Union[OptionKey, str]) -> ElementaryOptionValues:
        return self.get_value_for(key)

    def __len__(self) -> int:
        return len(self.options)

    def get_value_object_for(self, key: 'T.Union[OptionKey, str]') -> AnyOptionType:
        key = self.ensure_and_validate_key(key)
        potential = self.options.get(key, None)
        if self.is_project_option(key):
            assert key.subproject is not None
            if potential is None:
                raise KeyError(f'Tried to access nonexistant project option {key}.')
        else:
            if potential is None:
                parent_key = OptionKey(key.name, subproject=None, machine=key.machine)
                if parent_key not in self.options:
                    raise KeyError(f'Tried to access nonexistant project parent option {parent_key}.')
                # This is a global option but it can still have per-project
                # augment, so return the subproject key.
                return self.options[parent_key]
        return potential

    def get_value_object_and_value_for(self, key: OptionKey) -> T.Tuple[AnyOptionType, ElementaryOptionValues]:
        assert isinstance(key, OptionKey)
        key = self.ensure_and_validate_key(key)
        vobject = self.get_value_object_for(key)
        computed_value = vobject.value
        if key in self.augments:
            assert key.subproject is not None
            computed_value = self.augments[key]
        elif vobject.yielding:
            computed_value = vobject.parent.value
        return (vobject, computed_value)

    def option_has_value(self, key: OptionKey, value: ElementaryOptionValues) -> bool:
        vobject, current_value = self.get_value_object_and_value_for(key)
        return vobject.validate_value(value) == current_value

    def get_value_for(self, name: 'T.Union[OptionKey, str]', subproject: T.Optional[str] = None) -> ElementaryOptionValues:
        if isinstance(name, str):
            key = OptionKey(name, subproject)
        else:
            assert subproject is None
            key = name
        vobject, resolved_value = self.get_value_object_and_value_for(key)
        return resolved_value

    def add_system_option(self, key: T.Union[OptionKey, str], valobj: AnyOptionType) -> None:
        key = self.ensure_and_validate_key(key)
        if '.' in key.name:
            raise MesonException(f'Internal error: non-module option has a period in its name {key.name}.')
        self.add_system_option_internal(key, valobj)

    def add_system_option_internal(self, key: OptionKey, valobj: AnyOptionType) -> None:
        assert isinstance(valobj, UserOption)
        if not isinstance(valobj.name, str):
            assert isinstance(valobj.name, str)
        if key in self.options:
            return

        pval = self.pending_options.pop(key, None)
        if key.subproject:
            proj_key = key.evolve(subproject=None)
            self.add_system_option_internal(proj_key, valobj)
            if pval is not None:
                self.augments[key] = pval
        else:
            self.options[key] = valobj
            if pval is not None:
                self.set_option(key, pval)

    def add_compiler_option(self, language: str, key: T.Union[OptionKey, str], valobj: AnyOptionType) -> None:
        key = self.ensure_and_validate_key(key)
        if not key.name.startswith(language + '_'):
            raise MesonException(f'Internal error: all compiler option names must start with language prefix. ({key.name} vs {language}_)')
        self.add_system_option(key, valobj)

    def add_project_option(self, key: T.Union[OptionKey, str], valobj: AnyOptionType) -> None:
        key = self.ensure_and_validate_key(key)
        assert key.subproject is not None
        if key in self.options:
            raise MesonException(f'Internal error: tried to add a project option {key} that already exists.')
        if valobj.yielding and key.subproject:
            parent_key = key.as_root()
            try:
                parent_option = self.options[parent_key]
                # If parent object has different type, do not yield.
                # This should probably be an error.
                if type(parent_option) is type(valobj):
                    valobj.parent = parent_option
            except KeyError:
                # Subproject is set to yield, but top level
                # project does not have an option of the same
                pass
        valobj.yielding = bool(valobj.parent)

        self.options[key] = valobj
        self.project_options.add(key)
        pval = self.pending_options.pop(key, None)
        if pval is not None:
            self.set_option(key, pval)

    def add_module_option(self, modulename: str, key: T.Union[OptionKey, str], valobj: AnyOptionType) -> None:
        key = self.ensure_and_validate_key(key)
        if key.name.startswith('build.'):
            raise MesonException('FATAL internal error: somebody goofed option handling.')
        if not key.name.startswith(modulename + '.'):
            raise MesonException('Internal error: module option name {key.name} does not start with module prefix {modulename}.')
        self.add_system_option_internal(key, valobj)
        self.module_options.add(key)

    def sanitize_prefix(self, prefix: str) -> str:
        prefix = os.path.expanduser(prefix)
        if not os.path.isabs(prefix):
            raise MesonException(f'prefix value {prefix!r} must be an absolute path')
        if prefix.endswith('/') or prefix.endswith('\\'):
            # On Windows we need to preserve the trailing slash if the
            # string is of type 'C:\' because 'C:' is not an absolute path.
            if len(prefix) == 3 and prefix[1] == ':':
                pass
            # If prefix is a single character, preserve it since it is
            # the root directory.
            elif len(prefix) == 1:
                pass
            else:
                prefix = prefix[:-1]
        return prefix

    def sanitize_dir_option_value(self, prefix: str, option: OptionKey, value: T.Any) -> T.Any:
        '''
        If the option is an installation directory option, the value is an
        absolute path and resides within prefix, return the value
        as a path relative to the prefix. Otherwise, return it as is.

        This way everyone can do f.ex, get_option('libdir') and usually get
        the library directory relative to prefix, even though it really
        should not be relied upon.
        '''
        try:
            value = pathlib.PurePath(value)
        except TypeError:
            return value
        if option.name.endswith('dir') and value.is_absolute() and \
           option not in BUILTIN_DIR_NOPREFIX_OPTIONS:
            try:
                # Try to relativize the path.
                value = value.relative_to(prefix)
            except ValueError:
                # Path is not relative, let’s keep it as is.
                pass
            if '..' in value.parts:
                raise MesonException(
                    f"The value of the '{option}' option is '{value}' but "
                    "directory options are not allowed to contain '..'.\n"
                    f"If you need a path outside of the {prefix!r} prefix, "
                    "please use an absolute path."
                )
        # .as_posix() keeps the posix-like file separators Meson uses.
        return value.as_posix()

    def set_option(self, key: OptionKey, new_value: ElementaryOptionValues, first_invocation: bool = False) -> bool:
        changed = False
        error_key = key
        if error_key.subproject == '':
            error_key = error_key.evolve(subproject=None)

        if key.name == 'prefix':
            assert isinstance(new_value, str), 'for mypy'
            new_value = self.sanitize_prefix(new_value)
        elif self.is_builtin_option(key):
            prefix = self.get_value_for('prefix')
            assert isinstance(prefix, str), 'for mypy'
            new_value = self.sanitize_dir_option_value(prefix, key, new_value)

        try:
            opt = self.get_value_object_for(key)
        except KeyError:
            raise MesonException(f'Unknown option: "{error_key}".')

        if opt.deprecated is True:
            mlog.deprecation(f'Option "{error_key}" is deprecated')
        elif isinstance(opt.deprecated, list):
            for v in opt.listify(new_value):
                if v in opt.deprecated:
                    mlog.deprecation(f'Option "{error_key}" value {v!r} is deprecated')
        elif isinstance(opt.deprecated, dict):
            def replace(v: str) -> str:
                assert isinstance(opt.deprecated, dict) # No, Mypy can not tell this from two lines above
                newvalue = opt.deprecated.get(v)
                if newvalue is not None:
                    mlog.deprecation(f'Option "{error_key}" value {v!r} is replaced by {newvalue!r}')
                    return newvalue
                return v
            valarr = [replace(v) for v in opt.listify(new_value)]
            new_value = ','.join(valarr)
        elif isinstance(opt.deprecated, str):
            mlog.deprecation(f'Option "{error_key}" is replaced by {opt.deprecated!r}')
            # Change both this aption and the new one pointed to.
            changed |= self.set_option(key.evolve(name=opt.deprecated), new_value, first_invocation)

        new_value = opt.validate_value(new_value)
        if key in self.options:
            old_value = opt.value
            opt.set_value(new_value)
            opt.yielding = False
        else:
            assert key.subproject is not None
            old_value = self.augments.get(key, opt.value)
            self.augments[key] = new_value

        changed |= old_value != new_value
        if opt.readonly and changed and not first_invocation:
            raise MesonException(f'Tried to modify read only option "{error_key}"')

        if key.name == 'prefix' and first_invocation and changed:
            assert isinstance(old_value, str), 'for mypy'
            assert isinstance(new_value, str), 'for mypy'
            self.reset_prefixed_options(old_value, new_value)

        if changed and key.name == 'buildtype' and new_value != 'custom':
            assert isinstance(new_value, str), 'for mypy'
            optimization, debug = self.DEFAULT_DEPENDENTS[new_value]
            dkey = key.evolve(name='debug')
            optkey = key.evolve(name='optimization')
            self.set_option(dkey, debug, first_invocation)
            self.set_option(optkey, optimization, first_invocation)

        return changed

    def set_user_option(self, o: OptionKey, new_value: ElementaryOptionValues, first_invocation: bool = False) -> bool:
        if not self.is_cross and o.is_for_build():
            return False

        # This is complicated by the fact that a string can have two meanings:
        #
        # default_options: 'foo=bar'
        #
        # can be either
        #
        # A) a system option in which case the subproject is None
        # B) a project option, in which case the subproject is ''
        #
        # The key parsing function can not handle the difference between the two
        # and defaults to A.
        if o in self.options:
            return self.set_option(o, new_value, first_invocation)

        # could also be an augment...
        global_option = o.evolve(subproject=None)
        if o.subproject is not None and global_option in self.options:
            return self.set_option(o, new_value, first_invocation)

        if self.accept_as_pending_option(o, first_invocation=first_invocation):
            old_value = self.pending_options.get(o, None)
            self.pending_options[o] = new_value
            return old_value is None or str(old_value) != new_value
        elif o.subproject is None:
            o = o.as_root()
            return self.set_option(o, new_value, first_invocation)
        else:
            raise MesonException(f'Unknown option: "{o}".')

    def set_from_configure_command(self, D_args: T.Dict[OptionKey, T.Optional[str]]) -> bool:
        dirty = False
        for key, valstr in D_args.items():
            if valstr is not None:
                dirty |= self.set_user_option(key, valstr)
                continue

            if key in self.augments:
                del self.augments[key]
                dirty = True
            else:
                if key not in self.options:
                    raise MesonException(f"Unknown option: {key}")

                # TODO: For project options, "dropping an augment" means going
                # back to the superproject's value.  However, it's confusing
                # that -U does not simply remove the option from the stored
                # cmd_line_options.  This may cause "meson setup --wipe" to
                # have surprising behavior.  For this to work, UserOption
                # should only store the default value and the option values
                # should be stored with their source (project(), subproject(),
                # machine file, command line).  This way the effective value
                # can be easily recomputed.
                opt = self.get_value_object(key)
                dirty |= not opt.yielding and bool(opt.parent)
                opt.yielding = bool(opt.parent)
        return dirty

    def reset_prefixed_options(self, old_prefix: str, new_prefix: str) -> None:
        for optkey, prefix_mapping in BUILTIN_DIR_NOPREFIX_OPTIONS.items():
            valobj = self.options[optkey]
            new_value = valobj.value
            if new_prefix not in prefix_mapping:
                new_value = BUILTIN_OPTIONS[optkey].default
            else:
                if old_prefix in prefix_mapping:
                    # Only reset the value if it has not been changed from the default.
                    if prefix_mapping[old_prefix] == valobj.value:
                        new_value = prefix_mapping[new_prefix]
                else:
                    new_value = prefix_mapping[new_prefix]
            valobj.set_value(new_value)

    # FIXME, this should be removed.or renamed to "change_type_of_existing_object" or something like that
    def set_value_object(self, key: T.Union[OptionKey, str], new_object: AnyOptionType) -> None:
        key = self.ensure_and_validate_key(key)
        self.options[key] = new_object

    def get_value_object(self, key: T.Union[OptionKey, str]) -> AnyOptionType:
        key = self.ensure_and_validate_key(key)
        return self.options[key]

    def get_default_for_b_option(self, key: OptionKey) -> ElementaryOptionValues:
        assert self.is_base_option(key)
        try:
            return T.cast('ElementaryOptionValues', COMPILER_BASE_OPTIONS[key.evolve(subproject=None)].default)
        except KeyError:
            raise MesonBugException(f'Requested base option {key} which does not exist.')

    def remove(self, key: OptionKey) -> None:
        del self.options[key]
        try:
            self.project_options.remove(key)
        except KeyError:
            pass

    def __contains__(self, key: T.Union[str, OptionKey]) -> bool:
        key = self.ensure_and_validate_key(key)
        return key in self.options

    def __repr__(self) -> str:
        return repr(self.options)

    def keys(self) -> T.KeysView[OptionKey]:
        return self.options.keys()

    def values(self) -> T.ValuesView[AnyOptionType]:
        return self.options.values()

    def items(self) -> T.ItemsView['OptionKey', 'AnyOptionType']:
        return self.options.items()

    # FIXME: this method must be deleted and users moved to use "add_xxx_option"s instead.
    def update(self, **kwargs: AnyOptionType) -> None:
        self.options.update(**kwargs)

    def setdefault(self, k: OptionKey, o: AnyOptionType) -> AnyOptionType:
        return self.options.setdefault(k, o)

    def get(self, o: OptionKey, default: T.Optional[AnyOptionType] = None, **kwargs: T.Any) -> T.Optional[AnyOptionType]:
        return self.options.get(o, default, **kwargs)

    def is_project_option(self, key: OptionKey) -> bool:
        """Convenience method to check if this is a project option."""
        return key in self.project_options

    def is_per_machine_option(self, optname: OptionKey) -> bool:
        if optname.evolve(subproject=None, machine=MachineChoice.HOST) in BUILTIN_OPTIONS_PER_MACHINE:
            return True
        return self.is_compiler_option(optname)

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
        if isinstance(key, str):
            name: str = key
        else:
            name = key.name
        return name.startswith('backend_')

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

    def prefix_split_options(self, coll: OptionDict) -> T.Tuple[T.Optional[str], OptionDict]:
        prefix = None
        others_d: OptionDict = {}
        for k, v in coll.items():
            if k.name == 'prefix':
                if not isinstance(v, str):
                    raise MesonException('Incorrect type for prefix option (expected string)')
                prefix = v
            else:
                others_d[k] = v
        return (prefix, others_d)

    def first_handle_prefix(self,
                            project_default_options: OptionDict,
                            cmd_line_options: OptionDict,
                            machine_file_options: OptionDict) \
            -> T.Tuple[OptionDict, OptionDict, OptionDict]:
        # Copy to avoid later mutation
        nopref_machine_file_options = copy.copy(machine_file_options)

        prefix = None
        (possible_prefix, nopref_project_default_options) = self.prefix_split_options(project_default_options)
        prefix = prefix if possible_prefix is None else possible_prefix

        possible_prefixv = nopref_machine_file_options.pop(OptionKey('prefix'), None)
        assert possible_prefixv is None or isinstance(possible_prefixv, str), 'mypy: prefix from machine file was not a string?'
        prefix = prefix if possible_prefixv is None else possible_prefixv

        (possible_prefix, nopref_cmd_line_options) = self.prefix_split_options(cmd_line_options)
        prefix = prefix if possible_prefix is None else possible_prefix

        if prefix is not None:
            self.hard_reset_from_prefix(prefix)
        return (nopref_project_default_options, nopref_cmd_line_options, nopref_machine_file_options)

    def hard_reset_from_prefix(self, prefix: str) -> None:
        prefix = self.sanitize_prefix(prefix)
        for optkey, prefix_mapping in BUILTIN_DIR_NOPREFIX_OPTIONS.items():
            valobj = self.options[optkey]
            if prefix in prefix_mapping:
                new_value = prefix_mapping[prefix]
            else:
                _v = BUILTIN_OPTIONS[optkey].default
                assert isinstance(_v, str), 'for mypy'
                new_value = _v
            valobj.set_value(new_value)
        self.options[OptionKey('prefix')].set_value(prefix)

    def initialize_from_top_level_project_call(self,
                                               project_default_options_in: OptionDict,
                                               cmd_line_options_in: OptionDict,
                                               machine_file_options_in: OptionDict) -> None:
        (project_default_options, cmd_line_options, machine_file_options) = self.first_handle_prefix(project_default_options_in,
                                                                                                     cmd_line_options_in,
                                                                                                     machine_file_options_in)
        for key, valstr in project_default_options.items():
            # Due to backwards compatibility we ignore build-machine options
            # when building natively.
            if not self.is_cross and key.is_for_build():
                continue
            if key.subproject:
                # Subproject options from toplevel project() have low priority
                # and will be processed when the subproject is found
                self.pending_subproject_options[key] = valstr
            else:
                # Setting a project option with default_options
                # should arguably be a hard error; the default
                # value of project option should be set in the option
                # file, not in the project call.
                self.set_user_option(key, valstr, True)

        # ignore subprojects for now for machine file and command line
        # options; they are applied later
        for key, valstr in machine_file_options.items():
            # Due to backwards compatibility we ignore all build-machine options
            # when building natively.
            if not self.is_cross and key.is_for_build():
                continue
            if not key.subproject:
                self.set_user_option(key, valstr, True)
        for key, valstr in cmd_line_options.items():
            # Due to backwards compatibility we ignore all build-machine options
            # when building natively.
            if not self.is_cross and key.is_for_build():
                continue
            if not key.subproject:
                self.set_user_option(key, valstr, True)

    def accept_as_pending_option(self, key: OptionKey, first_invocation: bool = False) -> bool:
        # Some base options (sanitizers etc) might get added later.
        # Permitting them all is not strictly correct.
        if self.is_compiler_option(key):
            return True
        if first_invocation and self.is_backend_option(key):
            return True
        return (self.is_base_option(key) and
                key.evolve(subproject=None, machine=MachineChoice.HOST) in COMPILER_BASE_OPTIONS)

    def validate_cmd_line_options(self, cmd_line_options: OptionDict) -> None:
        unknown_options = []
        for key, valstr in cmd_line_options.items():
            if key in self.pending_options and not self.accept_as_pending_option(key):
                unknown_options.append(f'"{key}"')

        if unknown_options:
            keys = ', '.join(unknown_options)
            raise MesonException(f'Unknown options: {keys}')

    def initialize_from_subproject_call(self,
                                        subproject: str,
                                        spcall_default_options: OptionDict,
                                        project_default_options: OptionDict,
                                        cmd_line_options: OptionDict,
                                        machine_file_options: OptionDict) -> None:

        options: OptionDict = {}

        # project() default_options
        for key, valstr in project_default_options.items():
            if key.subproject == subproject:
                without_subp = key.evolve(subproject=None)
                raise MesonException(f'subproject name not needed in default_options; use "{without_subp}" instead of "{key}"')

            if key.subproject is None:
                key = key.evolve(subproject=subproject)
            options[key] = valstr

        # then global settings from machine file and command line
        # **but not if they are toplevel project options**
        for key, valstr in itertools.chain(machine_file_options.items(), cmd_line_options.items()):
            if key.subproject is None and not self.is_project_option(key.as_root()):
                subp_key = key.evolve(subproject=subproject)
                # just leave in place the value that was set for the toplevel project
                options.pop(subp_key, None)

        # augments from the toplevel project() default_options
        for key, valstr in self.pending_subproject_options.items():
            if key.subproject == subproject:
                options[key] = valstr

        # subproject() default_options
        for key, valstr in spcall_default_options.items():
            if key.subproject == subproject:
                without_subp = key.evolve(subproject=None)
                raise MesonException(f'subproject name not needed in default_options; use "{without_subp}" instead of "{key}"')

            if key.subproject is None:
                key = key.evolve(subproject=subproject)
            options[key] = valstr

        # then finally per project augments from machine file and command line
        for key, valstr in itertools.chain(machine_file_options.items(), cmd_line_options.items()):
            if key.subproject == subproject:
                options[key] = valstr

        # merge everything that has been computed above, while giving self.augments priority
        for key, valstr in options.items():
            if key.subproject != subproject:
                if key.subproject in self.subprojects and not self.option_has_value(key, valstr):
                    mlog.warning('option {key} is set in subproject {subproject} but has already been processed')
                    continue

                # Subproject options from project() will be processed when the subproject is found
                self.pending_subproject_options[key] = valstr
                continue

            self.pending_subproject_options.pop(key, None)
            self.pending_options.pop(key, None)
            if key not in self.augments:
                self.set_user_option(key, valstr, True)

        self.subprojects.add(subproject)

    def update_project_options(self, project_options: MutableKeyedOptionDictType, subproject: SubProject) -> None:
        for key, value in project_options.items():
            if key not in self.options:
                self.add_project_option(key, value)
                continue
            if key.subproject != subproject:
                raise MesonBugException(f'Tried to set an option for subproject {key.subproject} from {subproject}!')

            oldval = self.get_value_object(key)
            if type(oldval) is not type(value):
                self.set_option(key, value.value)
            elif choices_are_different(oldval, value):
                # If the choices have changed, use the new value, but attempt
                # to keep the old options. If they are not valid keep the new
                # defaults but warn.
                self.set_value_object(key, value)
                try:
                    value.set_value(oldval.value)
                except MesonException:
                    mlog.warning(f'Old value(s) of {key} are no longer valid, resetting to default ({value.value}).',
                                 fatal=False)

        # Find any extranious keys for this project and remove them
        potential_removed_keys = self.options.keys() - project_options.keys()
        for key in potential_removed_keys:
            if self.is_project_option(key) and key.subproject == subproject:
                self.remove(key)
