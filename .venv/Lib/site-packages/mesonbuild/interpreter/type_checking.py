# SPDX-License-Identifier: Apache-2.0
# Copyright © 2021-2025 Intel Corporation

"""Helpers for strict type checking."""

from __future__ import annotations
import itertools, os, re
import typing as T

from .. import compilers
from ..build import (CustomTarget, BuildTarget,
                     CustomTargetIndex, ExtractedObjects, GeneratedList, IncludeDirs,
                     BothLibraries, SharedLibrary, StaticLibrary, Jar, Executable, StructuredSources)
from ..options import OptionKey, UserFeatureOption
from ..dependencies import Dependency, DependencyMethods, InternalDependency
from ..interpreterbase.decorators import KwargInfo, ContainerTypeInfo, FeatureBroken, FeatureDeprecated
from ..mesonlib import (File, FileMode, MachineChoice, has_path_sep, listify, stringlistify,
                        EnvironmentVariables)
from ..programs import Program, ExternalProgram

# Helper definition for type checks that are `Optional[T]`
NoneType: T.Type[None] = type(None)

if T.TYPE_CHECKING:
    from typing_extensions import Literal

    from ..build import ObjectTypes, GeneratedTypes, BuildTargetTypes
    from ..interpreterbase import TYPE_var
    from ..options import ElementaryOptionValues
    from ..mesonlib import EnvInitValueType
    from ..interpreterbase.decorators import FeatureCheckBase

    _FullEnvInitValueType = T.Union[EnvironmentVariables, T.List[str], T.List[T.List[str]], EnvInitValueType, str, None]
    PkgConfigDefineType = T.Optional[T.Tuple[T.Tuple[str, str], ...]]
    SourcesVarargsType = T.List[T.Union[str, File, GeneratedTypes, StructuredSources, ExtractedObjects, BuildTarget]]


def in_set_validator(choices: T.Set[str]) -> T.Callable[[str], T.Optional[str]]:
    """Check that the choice given was one of the given set."""

    def inner(check: str) -> T.Optional[str]:
        if check not in choices:
            return f"must be one of {', '.join(sorted(choices))}, not {check}"
        return None

    return inner


def _language_validator(l: T.List[str]) -> T.Optional[str]:
    """Validate language keyword argument.

    Particularly for functions like `add_compiler()`, and `add_*_args()`
    """
    diff = {a.lower() for a in l}.difference(compilers.all_languages)
    if diff:
        return f'unknown languages: {", ".join(diff)}'
    return None


def _install_mode_validator(mode: T.List[T.Union[str, bool, int]]) -> T.Optional[str]:
    """Validate the `install_mode` keyword argument.

    This is a rather odd thing, it's a scalar, or an array of 3 values in the form:
    [(str | False), (str | int | False) = False, (str | int | False) = False]
    where the second and third components are not required and default to False.
    """
    if not mode:
        return None
    if True in mode:
        return 'components can only be permission strings, numbers, or False'
    if len(mode) > 3:
        return 'may have at most 3 elements'

    perms = mode[0]
    if not isinstance(perms, (str, bool)):
        return 'first component must be a permissions string or False'

    if isinstance(perms, str):
        if not len(perms) == 9:
            return ('permissions string must be exactly 9 characters in the form rwxr-xr-x,'
                    f' got {len(perms)}')
        for i in [0, 3, 6]:
            if perms[i] not in {'-', 'r'}:
                return f'permissions character {i+1} must be "-" or "r", not {perms[i]}'
        for i in [1, 4, 7]:
            if perms[i] not in {'-', 'w'}:
                return f'permissions character {i+1} must be "-" or "w", not {perms[i]}'
        for i in [2, 5]:
            if perms[i] not in {'-', 'x', 's', 'S'}:
                return f'permissions character {i+1} must be "-", "s", "S", or "x", not {perms[i]}'
        if perms[8] not in {'-', 'x', 't', 'T'}:
            return f'permission character 9 must be "-", "t", "T", or "x", not {perms[8]}'

        if len(mode) >= 2 and not isinstance(mode[1], (int, str, bool)):
            return 'second component can only be a string, number, or False'
        if len(mode) >= 3 and not isinstance(mode[2], (int, str, bool)):
            return 'third component can only be a string, number, or False'

    return None


def _install_mode_convertor(mode: T.Optional[T.List[T.Union[str, bool, int]]]) -> FileMode:
    """Convert the DSL form of the `install_mode` keyword argument to `FileMode`"""

    if not mode:
        return FileMode()

    # This has already been validated by the validator. False denotes "use
    # default". mypy is totally incapable of understanding it, because
    # generators clobber types via homogeneous return. But also we *must*
    # convert the first element different from the rest
    m1 = mode[0] if isinstance(mode[0], str) else None
    rest = (m if isinstance(m, (str, int)) else None for m in mode[1:])

    return FileMode(m1, *rest)


def _lower_strlist(input: T.List[str]) -> T.List[str]:
    """Lower a list of strings.

    mypy (but not pyright) gets confused about using a lambda as the convertor function
    """
    return [i.lower() for i in input]


def _validate_shlib_version(val: T.Optional[str]) -> T.Optional[str]:
    if val is not None and not re.fullmatch(r'[0-9]+(\.[0-9]+){0,2}', val):
        return (f'Invalid Shared library version "{val}". '
                'Must be of the form X.Y.Z where all three are numbers. Y and Z are optional.')
    return None


def variables_validator(contents: T.Union[str, T.List[str], T.Dict[str, str]]) -> T.Optional[str]:
    if isinstance(contents, str):
        contents = [contents]
    if isinstance(contents, dict):
        variables = contents
    else:
        variables = {}
        for v in contents:
            try:
                key, val = v.split('=', 1)
            except ValueError:
                return f'variable {v!r} must have a value separated by equals sign.'
            variables[key.strip()] = val.strip()
    for k, v in variables.items():
        if not k:
            return 'empty variable name'
        if any(c.isspace() for c in k):
            return f'invalid whitespace in variable name {k!r}'
    return None


def variables_convertor(contents: T.Union[str, T.List[str], T.Dict[str, str]]) -> T.Dict[str, str]:
    if isinstance(contents, str):
        contents = [contents]
    if isinstance(contents, dict):
        return contents
    variables = {}
    for v in contents:
        key, val = v.split('=', 1)
        variables[key.strip()] = val.strip()
    return variables


NATIVE_KW = KwargInfo(
    'native', bool,
    default=False,
    convertor=lambda n: MachineChoice.BUILD if n else MachineChoice.HOST)

LANGUAGE_KW = KwargInfo(
    'language', ContainerTypeInfo(list, str, allow_empty=False),
    listify=True,
    required=True,
    validator=_language_validator,
    convertor=_lower_strlist)

INSTALL_MODE_KW: KwargInfo[T.List[T.Union[str, bool, int]]] = KwargInfo(
    'install_mode',
    ContainerTypeInfo(list, (str, bool, int)),
    listify=True,
    default=[],
    validator=_install_mode_validator,
    convertor=_install_mode_convertor,
)

REQUIRED_KW: KwargInfo[T.Union[bool, UserFeatureOption]] = KwargInfo(
    'required',
    (bool, UserFeatureOption),
    default=True,
    # TODO: extract_required_kwarg could be converted to a convertor
)

DISABLER_KW: KwargInfo[bool] = KwargInfo('disabler', bool, default=False)

def _env_validator(value: T.Union[EnvironmentVariables, T.List['TYPE_var'], T.Dict[str, 'TYPE_var'], str, None],
                   only_dict_str: bool = True) -> T.Optional[str]:
    def _splitter(v: str) -> T.Optional[str]:
        split = v.split('=', 1)
        if len(split) == 1:
            return f'"{v}" is not two string values separated by an "="'
        return None

    if isinstance(value, str):
        v = _splitter(value)
        if v is not None:
            return v
    elif isinstance(value, list):
        for i in listify(value):
            if not isinstance(i, str):
                return f"All array elements must be a string, not {i!r}"
            v = _splitter(i)
            if v is not None:
                return v
    elif isinstance(value, dict):
        # We don't need to spilt here, just do the type checking
        for k, dv in value.items():
            if only_dict_str:
                if any(i for i in listify(dv) if not isinstance(i, str)):
                    return f"Dictionary element {k} must be a string or list of strings not {dv!r}"
            elif isinstance(dv, list):
                if any(not isinstance(i, str) for i in dv):
                    return f"Dictionary element {k} must be a string, bool, integer or list of strings, not {dv!r}"
            elif not isinstance(dv, (str, bool, int)):
                return f"Dictionary element {k} must be a string, bool, integer or list of strings, not {dv!r}"
    # We know that otherwise we have an EnvironmentVariables object or None, and
    # we're okay at this point
    return None

def _options_validator(value: T.Union[EnvironmentVariables, T.List['TYPE_var'], T.Dict[str, 'TYPE_var'], str, None]) -> T.Optional[str]:
    # Reusing the env validator is a little overkill, but nicer than duplicating the code
    return _env_validator(value, only_dict_str=False)

def split_equal_string(input: str) -> T.Tuple[str, str]:
    """Split a string in the form `x=y`

    This assumes that the string has already been validated to split properly.
    """
    a, b = input.split('=', 1)
    return (a, b)

# Split _env_convertor() and env_convertor_with_method() to make mypy happy.
# It does not want extra arguments in KwargInfo convertor callable.
def env_convertor_with_method(value: _FullEnvInitValueType,
                              init_method: Literal['set', 'prepend', 'append'] = 'set',
                              separator: str = os.pathsep) -> EnvironmentVariables:
    if isinstance(value, str):
        return EnvironmentVariables(dict([split_equal_string(value)]), init_method, separator)
    elif isinstance(value, list):
        return EnvironmentVariables(dict(split_equal_string(v) for v in listify(value)), init_method, separator)
    elif isinstance(value, dict):
        return EnvironmentVariables({k: listify(dv) for k, dv in value.items()}, init_method, separator)
    elif value is None:
        return EnvironmentVariables()
    return value

def _env_convertor(value: _FullEnvInitValueType) -> EnvironmentVariables:
    return env_convertor_with_method(value)

ENV_KW: KwargInfo[T.Union[EnvironmentVariables, T.List, T.Dict, str, None]] = KwargInfo(
    'env',
    (EnvironmentVariables, list, dict, str, NoneType),
    validator=_env_validator,
    convertor=_env_convertor,
)

DEPFILE_KW: KwargInfo[T.Optional[str]] = KwargInfo(
    'depfile',
    (str, type(None)),
    validator=lambda x: 'Depfile must be a plain filename with a subdirectory' if has_path_sep(x) else None
)

DEPENDS_KW: KwargInfo[T.List[BuildTargetTypes]] = KwargInfo(
    'depends',
    ContainerTypeInfo(list, (BuildTarget, CustomTarget, CustomTargetIndex)),
    listify=True,
    default=[],
    since_values={CustomTargetIndex: '1.5.0'},
)

DEPEND_FILES_KW: KwargInfo[T.List[T.Union[str, File]]] = KwargInfo(
    'depend_files',
    ContainerTypeInfo(list, (File, str)),
    listify=True,
    default=[],
)

COMMAND_KW: KwargInfo[T.List[T.Union[str, BuildTargetTypes, Program, File]]] = KwargInfo(
    'command',
    ContainerTypeInfo(list, (str, BuildTarget, CustomTarget, CustomTargetIndex, Program, File), allow_empty=False),
    required=True,
    listify=True,
    default=[],
)


def _override_options_convertor(raw: T.Union[str, T.List[str], T.Dict[str, ElementaryOptionValues]]) -> T.Dict[str, ElementaryOptionValues]:
    if isinstance(raw, dict):
        return raw
    raw = stringlistify(raw)
    output: T.Dict[str, ElementaryOptionValues] = {}
    for each in raw:
        k, v = split_equal_string(each)
        output[k] = v
    return output

OVERRIDE_OPTIONS_KW: KwargInfo[T.Union[str, T.List[str], T.Dict[str, ElementaryOptionValues]]] = KwargInfo(
    'override_options',
    (str, ContainerTypeInfo(list, str), ContainerTypeInfo(dict, (str, int, bool, list))),
    default={},
    validator=_options_validator,
    convertor=_override_options_convertor,
    since_values={dict: '1.2.0'},
)


def _output_validator(outputs: T.List[str]) -> T.Optional[str]:
    output_set = set(outputs)
    if len(output_set) != len(outputs):
        seen = set()
        for el in outputs:
            if el in seen:
                return f"contains {el!r} multiple times, but no duplicates are allowed."
            seen.add(el)
    for i in outputs:
        if i == '':
            return 'Output must not be empty.'
        elif i.strip() == '':
            return 'Output must not consist only of whitespace.'
        elif has_path_sep(i):
            return f'Output {i!r} must not contain a path segment.'
        elif '@INPUT' in i:
            return f'output {i!r} contains "@INPUT", which is invalid. Did you mean "@PLAINNAME@" or "@BASENAME@?'

    return None

MULTI_OUTPUT_KW: KwargInfo[T.List[str]] = KwargInfo(
    'output',
    ContainerTypeInfo(list, str, allow_empty=False),
    listify=True,
    required=True,
    default=[],
    validator=_output_validator,
)

OUTPUT_KW: KwargInfo[str] = KwargInfo(
    'output',
    str,
    required=True,
    validator=lambda x: _output_validator([x])
)

CT_INPUT_KW: KwargInfo[T.List[T.Union[str, File, ExternalProgram, BuildTarget, GeneratedTypes, ExtractedObjects]]] = KwargInfo(
    'input',
    ContainerTypeInfo(list, (str, File, ExternalProgram, BuildTarget, CustomTarget, CustomTargetIndex, ExtractedObjects, GeneratedList)),
    listify=True,
    default=[],
)

CT_INSTALL_TAG_KW: KwargInfo[T.List[T.Union[str, bool]]] = KwargInfo(
    'install_tag',
    ContainerTypeInfo(list, (str, bool)),
    listify=True,
    default=[],
    since='0.60.0',
    convertor=lambda x: [y if isinstance(y, str) else None for y in x],
)

INSTALL_TAG_KW: KwargInfo[T.Optional[str]] = KwargInfo('install_tag', (str, NoneType))

INSTALL_FOLLOW_SYMLINKS: KwargInfo[T.Optional[bool]] = KwargInfo(
    'follow_symlinks',
    (bool, NoneType),
    since='1.3.0',
)

INSTALL_KW = KwargInfo('install', bool, default=False)

CT_INSTALL_DIR_KW: KwargInfo[T.List[T.Union[str, Literal[False]]]] = KwargInfo(
    'install_dir',
    ContainerTypeInfo(list, (str, bool)),
    listify=True,
    default=[],
    validator=lambda x: 'must be `false` if boolean' if True in x else None,
)

CT_BUILD_BY_DEFAULT: KwargInfo[T.Optional[bool]] = KwargInfo('build_by_default', (bool, type(None)), since='0.40.0')

CT_BUILD_ALWAYS: KwargInfo[T.Optional[bool]] = KwargInfo(
    'build_always', (bool, NoneType),
    deprecated='0.47.0',
    deprecated_message='combine build_by_default and build_always_stale instead.',
)

CT_BUILD_ALWAYS_STALE: KwargInfo[T.Optional[bool]] = KwargInfo(
    'build_always_stale', (bool, NoneType),
    since='0.47.0',
)

INSTALL_DIR_KW: KwargInfo[T.Optional[str]] = KwargInfo('install_dir', (str, NoneType))

INCLUDE_DIRECTORIES: KwargInfo[T.List[T.Union[str, IncludeDirs]]] = KwargInfo(
    'include_directories',
    ContainerTypeInfo(list, (str, IncludeDirs)),
    listify=True,
    default=[],
)

def _default_options_convertor(raw: T.Union[str, T.List[str], T.Dict[str, ElementaryOptionValues]]) -> T.Dict[OptionKey, ElementaryOptionValues]:
    d = _override_options_convertor(raw)
    return {OptionKey.from_string(k): v for k, v in d.items()}

DEFAULT_OPTIONS = OVERRIDE_OPTIONS_KW.evolve(
        name='default_options',
        convertor=_default_options_convertor)

ENV_METHOD_KW = KwargInfo('method', str, default='set', since='0.62.0',
                          validator=in_set_validator({'set', 'prepend', 'append'}))

ENV_SEPARATOR_KW = KwargInfo('separator', str, default=os.pathsep)

DEPENDENCIES_KW: KwargInfo[T.List[Dependency]] = KwargInfo(
    'dependencies',
    # InternalDependency is a subclass of Dependency, but we want to
    # print it in error messages
    ContainerTypeInfo(list, (Dependency, InternalDependency)),
    listify=True,
    default=[],
    extra_types={
        BuildTarget: lambda arg: f'Tried to use a build_target "{T.cast("BuildTarget", arg).name}" as a dependency. This should be in `link_with` or `link_whole` instead.',
    },
)

D_MODULE_VERSIONS_KW: KwargInfo[T.List[T.Union[str, int]]] = KwargInfo(
    'd_module_versions',
    ContainerTypeInfo(list, (str, int)),
    listify=True,
    default=[],
)

_LINK_WITH_ERROR = 'Dependency and external_library objects must go in the "dependencies" keyword argument'

def _link_with_validator(values: T.List[T.Union[BothLibraries, SharedLibrary, StaticLibrary,
                                                CustomTarget, CustomTargetIndex, Jar, Executable,
                                                ]]
                         ) -> T.Optional[str]:
    for value in values:
        if not value.is_linkable_target():
            return f'Link target "{value!s}" is not linkable'
    return None

# Allow Dependency for the better error message? But then in other cases it will list this as one of the allowed types!
LINK_WITH_KW: KwargInfo[T.List[T.Union[BothLibraries, SharedLibrary, StaticLibrary, CustomTarget, CustomTargetIndex, Jar, Executable]]] = KwargInfo(
    'link_with',
    ContainerTypeInfo(list, (BothLibraries, SharedLibrary, StaticLibrary, CustomTarget, CustomTargetIndex, Jar, Executable)),
    listify=True,
    default=[],
    extra_types={Dependency: lambda _: _LINK_WITH_ERROR},
    validator=_link_with_validator,
)

def link_whole_validator(values: T.List[T.Union[StaticLibrary, CustomTarget, CustomTargetIndex]]) -> T.Optional[str]:
    for l in values:
        if isinstance(l, (CustomTarget, CustomTargetIndex)) and l.links_dynamically():
            return f'{type(l).__name__} returning a shared library is not allowed'
        if not l.is_linkable_target():
            return f'Link target "{l!s}" is not linkable'
    return None

LINK_WHOLE_KW: KwargInfo[T.List[T.Union[BothLibraries, StaticLibrary, CustomTarget, CustomTargetIndex]]] = KwargInfo(
    'link_whole',
    ContainerTypeInfo(list, (BothLibraries, StaticLibrary, CustomTarget, CustomTargetIndex)),
    listify=True,
    default=[],
    validator=link_whole_validator,
    extra_types={Dependency: lambda _: _LINK_WITH_ERROR}
)

DEPENDENCY_SOURCES_KW: KwargInfo[T.List[T.Union[str, File, GeneratedTypes]]] = KwargInfo(
    'sources',
    ContainerTypeInfo(list, (str, File, CustomTarget, CustomTargetIndex, GeneratedList)),
    listify=True,
    default=[],
)

SOURCES_VARARGS = (str, File, CustomTarget, CustomTargetIndex, GeneratedList, StructuredSources, ExtractedObjects, BuildTarget)

BT_SOURCES_KW: KwargInfo[SourcesVarargsType] = KwargInfo(
    'sources',
    (NoneType, ContainerTypeInfo(list, SOURCES_VARARGS)),
    listify=True,
    default=[],
)

VARIABLES_KW: KwargInfo[T.Dict[str, str]] = KwargInfo(
    'variables',
    # str is listified by validator/convertor, cannot use listify=True here because
    # that would listify dict too.
    (str, ContainerTypeInfo(list, str), ContainerTypeInfo(dict, str)), # type: ignore
    validator=variables_validator,
    convertor=variables_convertor,
    default={},
)

PRESERVE_PATH_KW: KwargInfo[bool] = KwargInfo('preserve_path', bool, default=False, since='0.63.0')

def suite_convertor(suite: T.List[str]) -> T.List[str]:
    # Ensure we always have at least one suite.
    if not suite:
        return ['']
    return suite

TEST_KWS_NO_ARGS: T.List[KwargInfo] = [
    KwargInfo('should_fail', (bool, NoneType), deprecated='1.11.0', deprecated_message='Use expected_fail instead of should_fail'),
    KwargInfo('expected_fail', (bool, NoneType), since='1.11.0'),
    KwargInfo('expected_exitcode', (int, NoneType), since='1.11.0'),
    KwargInfo('timeout', int, default=30),
    KwargInfo('workdir', (str, NoneType), default=None,
              validator=lambda x: 'must be an absolute path' if not os.path.isabs(x) else None),
    KwargInfo('protocol', str,
              default='exitcode',
              validator=in_set_validator({'exitcode', 'tap', 'gtest', 'rust'}),
              since_values={'gtest': '0.55.0', 'rust': '0.57.0'}),
    KwargInfo('priority', int, default=0, since='0.52.0'),
    # TODO: env needs reworks of the way the environment variable holder itself works probably
    ENV_KW,
    DEPENDS_KW.evolve(since='0.46.0'),
    KwargInfo('suite', ContainerTypeInfo(list, str), listify=True, default=[], convertor=suite_convertor),
    KwargInfo('verbose', bool, default=False, since='0.62.0'),
]

TEST_KWS: T.List[KwargInfo] = TEST_KWS_NO_ARGS + [
    KwargInfo('args', ContainerTypeInfo(list, (str, File, BuildTarget, CustomTarget, CustomTargetIndex, Program)),
              listify=True, default=[]),
]

# Cannot have a default value because we need to check that rust_crate_type and
# rust_abi are mutually exclusive.
RUST_CRATE_TYPE_KW: KwargInfo[T.Union[str, None]] = KwargInfo(
    'rust_crate_type', (str, NoneType),
    since='0.42.0',
    since_values={'proc-macro': '0.62.0'},
    deprecated='1.3.0',
    deprecated_message='Use rust_abi or rust.proc_macro() instead.',
    validator=in_set_validator({'bin', 'lib', 'rlib', 'dylib', 'cdylib', 'staticlib', 'proc-macro'}))

RUST_ABI_KW: KwargInfo[T.Union[str, None]] = KwargInfo(
    'rust_abi', (str, NoneType),
    since='1.3.0',
    validator=in_set_validator({'rust', 'c'}))

_VS_MODULE_DEFS_KW: KwargInfo[T.Optional[T.Union[str, File, CustomTarget, CustomTargetIndex]]] = KwargInfo(
    'vs_module_defs',
    (str, File, CustomTarget, CustomTargetIndex, NoneType),
    since_values={CustomTargetIndex: '1.3.0'}
)

_BASE_LANG_KW: KwargInfo[T.List[str]] = KwargInfo(
    'UNKNOWN',
    ContainerTypeInfo(list, (str)),
    listify=True,
    default=[],
)

_LANGUAGE_KWS: T.List[KwargInfo[T.List[str]]] = [
    _BASE_LANG_KW.evolve(name=f'{lang}_args')
    for lang in compilers.all_languages - {'rust', 'vala', 'java'}
]
# Cannot use _BASE_LANG_KW here because Vala is special for types
_LANGUAGE_KWS.append(KwargInfo(
    'vala_args', ContainerTypeInfo(list, (str, File)), listify=True, default=[]))
_LANGUAGE_KWS.append(_BASE_LANG_KW.evolve(name='rust_args', since='0.41.0'))

# We need this deprecated values more than the non-deprecated values. So we'll evolve them out elsewhere.
_JAVA_LANG_KW: KwargInfo[T.List[str]] = _BASE_LANG_KW.evolve(
    name='java_args',
    deprecated='1.3.0',
    deprecated_message='This does not, and never has, done anything. It should be removed'
)

def _objects_validator(vals: T.List[ObjectTypes]) -> T.Optional[str]:
    non_objects: T.List[str] = []

    for val in vals:
        if isinstance(val, (str, File, ExtractedObjects)):
            continue
        else:
            non_objects.extend(o for o in val.get_outputs() if not compilers.is_object(o))

    if non_objects:
        return f'{", ".join(non_objects)!r} are not objects'

    return None


def _target_install_feature_validator(val: object) -> T.Iterable[FeatureCheckBase]:
    # due to lack of type checking, these are "allowed" for legacy reasons
    if not isinstance(val, bool):
        yield FeatureBroken('install kwarg with non-boolean value', '1.3.0',
                            'This was never intended to work, and is essentially the same as using `install: true` regardless of value.')


def _target_install_convertor(val: object) -> bool:
    return bool(val)


def _extra_files_validator(args: T.List[T.Union[File, str]]) -> T.Optional[str]:
    generated = [a for a in args if isinstance(a, File) and a.is_built]
    if generated:
        return 'extra_files contains generated files: {}'.format(', '.join(f"{f.fname}" for f in generated))
    return None


def _bt_install_dir_deprecated(args: T.List[T.Union[str, bool]]) -> T.Iterator[FeatureCheckBase]:
    if len(args) > 1:
        yield FeatureDeprecated('passing more than one argument to install_dir', '1.11.0',
                                'use the install_vala_* arguments instead')


# Applies to all build_target like classes
_ALL_TARGET_KWS: T.List[KwargInfo] = [
    OVERRIDE_OPTIONS_KW,
    KwargInfo('build_by_default', bool, default=True, since='0.38.0'),
    DEPENDENCIES_KW,
    KwargInfo(
        'extra_files',
        ContainerTypeInfo(list, (str, File)),
        default=[],
        listify=True,
        validator=_extra_files_validator,
    ),
    INCLUDE_DIRECTORIES.evolve(since_values={ContainerTypeInfo(list, str): '0.50.0'}),
    KwargInfo(
        'install',
        object,
        default=False,
        convertor=_target_install_convertor,
        feature_validator=_target_install_feature_validator,
    ),
    INSTALL_MODE_KW,
    INSTALL_TAG_KW,
    KwargInfo(
        'install_dir',
        ContainerTypeInfo(list, (str, bool)),
        default=[],
        listify=True,
        feature_validator=_bt_install_dir_deprecated,
    ),
    KwargInfo('implicit_include_directories', bool, default=True, since='0.42.0'),
    LINK_WITH_KW.evolve(
        as_default=[('', ('1.11.0', "Replace an empty string with an empty array: `link_with : ''` -> `link_with : []`"))],
    ),
    NATIVE_KW,
    KwargInfo('resources', ContainerTypeInfo(list, str), default=[], listify=True),
    KwargInfo(
        'objects',
        ContainerTypeInfo(list, (str, File, CustomTarget, CustomTargetIndex, GeneratedList, ExtractedObjects)),
        listify=True,
        default=[],
        validator=_objects_validator,
        since_values={
            ContainerTypeInfo(list, (GeneratedList, CustomTarget, CustomTargetIndex)):
                ('1.1.0', 'generated sources as positional "objects" arguments')
        },
    ),
    KwargInfo('build_subdir', str, default='', since='1.10.0')
]


def _name_validator(arg: T.Optional[T.Union[str, T.List]]) -> T.Optional[str]:
    if isinstance(arg, list) and arg:
        return 'must be empty when passed as an array to signify the default value.'
    return None


def _name_suffix_validator(arg: T.Optional[T.Union[str, T.List]]) -> T.Optional[str]:
    if arg == '':
        return 'must not be a empty string. An empty array may be passed if you want Meson to use the default behavior.'
    return _name_validator(arg)


_NAME_PREFIX_KW: KwargInfo[T.Optional[T.Union[str, T.List]]] = KwargInfo(
    'name_prefix',
    (str, NoneType, list),
    validator=_name_validator,
    convertor=lambda x: None if isinstance(x, list) else x,
)


def _pch_validator(args: T.List[str]) -> T.Optional[str]:
    num_args = len(args)
    if num_args == 1:
        if not compilers.is_header(args[0]):
            return f'PCH argument {args[0]} is not a header.'
    elif num_args == 2:
        if compilers.is_header(args[0]):
            if not compilers.is_source(args[1]):
                return 'PCH definition must contain one header and at most one source.'
        elif compilers.is_source(args[0]):
            if not compilers.is_header(args[1]):
                return 'PCH definition must contain one header and at most one source.'
        else:
            return f'PCH argument {args[0]} has neither a known header or code extension.'

        if os.path.dirname(args[0]) != os.path.dirname(args[1]):
            return 'PCH files must be stored in the same folder.'
    elif num_args > 2:
        return 'A maximum of two elements are allowed for PCH arguments'
    if num_args >= 1 and not has_path_sep(args[0]):
        return f'PCH header {args[0]} must not be in the same directory as source files'
    if num_args == 2 and not has_path_sep(args[1]):
        return f'PCH source {args[0]} must not be in the same directory as source files'
    return None


def _pch_feature_validator(args: T.List[str]) -> T.Iterable[FeatureCheckBase]:
    if len(args) > 1:
        yield FeatureDeprecated('PCH source files', '0.50.0', 'Only a single header file should be used.')


def _pch_convertor(args: T.List[str]) -> T.Optional[T.Tuple[str, T.Optional[str]]]:
    num_args = len(args)

    if num_args == 1:
        return (args[0], None)

    if num_args == 2:
        if compilers.is_source(args[0]):
            # Flip so that we always have [header, src]
            return (args[1], args[0])
        return (args[0], args[1])

    return None


_PCH_ARGS: KwargInfo[T.List[str]] = KwargInfo(
    'pch',
    ContainerTypeInfo(list, str),
    listify=True,
    default=[],
    validator=_pch_validator,
    feature_validator=_pch_feature_validator,
    convertor=_pch_convertor,
)


LINK_ARGS_KW: KwargInfo[T.List[str]] = KwargInfo(
    'link_args',
    ContainerTypeInfo(list, str),
    default=[],
    listify=True,
    as_default=[('', ('1.10.1', "Replace an empty string with an empty array: `link_args : ''` -> `link_args : []`"))],
)


# Applies to all build_target classes except jar
_BUILD_TARGET_KWS: T.List[KwargInfo] = [
    *_ALL_TARGET_KWS,
    *_LANGUAGE_KWS,
    BT_SOURCES_KW,
    INCLUDE_DIRECTORIES.evolve(name='d_import_dirs'),
    LINK_ARGS_KW,
    LINK_WHOLE_KW.evolve(
        as_default=[('', ('1.11.0', "Replace an empty string with an empty array: `link_whole : ''` -> `link_whole : []`"))],
    ),
    _NAME_PREFIX_KW,
    _NAME_PREFIX_KW.evolve(name='name_suffix', validator=_name_suffix_validator),
    RUST_CRATE_TYPE_KW,
    _PCH_ARGS.evolve(name='c_pch'),
    _PCH_ARGS.evolve(name='cpp_pch'),
    KwargInfo('d_debug', ContainerTypeInfo(list, (str, int)), default=[], listify=True),
    D_MODULE_VERSIONS_KW,
    KwargInfo('d_unittest', bool, default=False),
    KwargInfo(
        'rust_dependency_map',
        ContainerTypeInfo(dict, str),
        default={},
        since='1.2.0',
    ),
    KwargInfo('swift_interoperability_mode', str, default='c', validator=in_set_validator({'c', 'cpp'}), since='1.9.0'),
    KwargInfo('swift_module_name', str, default='', since='1.9.0'),
    KwargInfo('build_rpath', str, default='', since='0.42.0'),
    KwargInfo(
        'gnu_symbol_visibility',
        str,
        default='',
        validator=in_set_validator({'', 'default', 'internal', 'hidden', 'protected', 'inlineshidden'}),
        since='0.48.0',
    ),
    KwargInfo('install_rpath', str, default=''),
    KwargInfo(
        'link_depends',
        ContainerTypeInfo(list, (str, File, CustomTarget, CustomTargetIndex, BuildTarget)),
        default=[],
        listify=True,
    ),
    KwargInfo(
        'link_language',
        (str, NoneType),
        # Yes, neither mypy no pyright can figure out that `set[literal[str]]``
        # is a subset of `set[str]`
        validator=in_set_validator(T.cast('T.Set[str]', compilers.all_languages)),
        since='0.51.0',
    ),
    KwargInfo('vala_gir', (str, NoneType)),
    KwargInfo('install_vala_gir', (str, bool, NoneType), since='1.11.0'),
    KwargInfo('vala_header', (str, NoneType)),
    KwargInfo('install_vala_header', (str, bool, NoneType), since='1.11.0'),
    KwargInfo('vala_vapi', (str, NoneType)),
    KwargInfo(
        'link_early_args',
        (ContainerTypeInfo(list, str), NoneType),
        listify=True,
        default=[],
        since='1.11.0',
    ),
    KwargInfo('install_vala_vapi', (str, bool, NoneType), since='1.11.0'),
]

def _validate_win_subsystem(value: T.Optional[str]) -> T.Optional[str]:
    if value is not None:
        if re.fullmatch(r'(boot_application|console|efi_application|efi_boot_service_driver|efi_rom|efi_runtime_driver|native|posix|windows)(,\d+(\.\d+)?)?', value) is None:
            return f'Invalid value for win_subsystem: {value}.'
    return None


def _validate_darwin_versions(darwin_versions: T.List[T.Union[str, int]]) -> T.Optional[str]:
    if len(darwin_versions) > 2:
        return f"Must contain between 0 and 2 elements, not {len(darwin_versions)}"
    if len(darwin_versions) == 1:
        darwin_versions = 2 * darwin_versions
    for v in darwin_versions:
        if isinstance(v, int):
            v = str(v)
        if not re.fullmatch(r'[0-9]+(\.[0-9]+){0,2}', v):
            return 'must be X.Y.Z where X, Y, Z are numbers, and Y and Z are optional'
        try:
            parts = v.split('.')
        except ValueError:
            return f'badly formed value: "{v}, not in X.Y.Z form'
        if len(parts) in {1, 2, 3} and int(parts[0]) > 65535:
            return 'must be X.Y.Z where X is [0, 65535] and Y, Z are optional'
        if len(parts) in {2, 3} and int(parts[1]) > 255:
            return 'must be X.Y.Z where Y is [0, 255] and Y, Z are optional'
        if len(parts) == 3 and int(parts[2]) > 255:
            return 'must be X.Y.Z where Z is [0, 255] and Y, Z are optional'
    return None


def _convert_darwin_versions(val: T.List[T.Union[str, int]]) -> T.Optional[T.Tuple[str, str]]:
    if not val:
        return None
    elif len(val) == 1:
        v = str(val[0])
        return (v, v)
    return (str(val[0]), str(val[1]))


_DARWIN_VERSIONS_KW: KwargInfo[T.List[T.Union[str, int]]] = KwargInfo(
    'darwin_versions',
    ContainerTypeInfo(list, (str, int)),
    default=[],
    listify=True,
    validator=_validate_darwin_versions,
    convertor=_convert_darwin_versions,
    since='0.48.0',
)

# Arguments exclusive to Executable. These are separated to make integrating
# them into build_target easier
EXCLUSIVE_EXECUTABLE_KWS: T.List[KwargInfo] = [
    KwargInfo('export_dynamic', (bool, NoneType), since='0.45.0'),
    KwargInfo('gui_app', (bool, NoneType), deprecated='0.56.0', deprecated_message="Use 'win_subsystem' instead"),
    KwargInfo(
        'implib',
        (bool, str, NoneType),
        since='0.42.0',
        deprecated_values={
            bool: ('1.10.0', 'Use "export_dynamic" keyword instead'),
        },
    ),
    KwargInfo('pie', (bool, NoneType)),
    KwargInfo(
        'win_subsystem',
        (str, NoneType),
        convertor=lambda x: x.lower() if isinstance(x, str) else None,
        validator=_validate_win_subsystem,
    ),
    KwargInfo(
        'android_exe_type',
        (str, NoneType),
        validator=in_set_validator({'application', 'executable'}),
        since='1.8.0'
    ),
]

# The total list of arguments used by Executable
EXECUTABLE_KWS = [
    *_BUILD_TARGET_KWS,
    *EXCLUSIVE_EXECUTABLE_KWS,
    _VS_MODULE_DEFS_KW.evolve(since='1.3.0', since_values=None),
    _JAVA_LANG_KW,
]

# Arguments exclusive to library types
_EXCLUSIVE_LIB_KWS: T.List[KwargInfo] = [
    RUST_ABI_KW,
]

# Arguments exclusive to StaticLibrary. These are separated to make integrating
# them into build_target easier
_EXCLUSIVE_STATIC_LIB_KWS: T.List[KwargInfo] = [
    KwargInfo('prelink', bool, default=False, since='0.57.0'),
    KwargInfo('pic', (bool, NoneType), since='0.36.0'),
]

# The total list of arguments used by StaticLibrary
STATIC_LIB_KWS = [
    *_BUILD_TARGET_KWS,
    *_EXCLUSIVE_STATIC_LIB_KWS,
    *_EXCLUSIVE_LIB_KWS,
    _JAVA_LANG_KW,
]

def _shortname_validator(shortname: T.Optional[str]) -> T.Optional[str]:
    if shortname is not None and len(shortname) > 8:
        return 'must have a maximum of 8 characters'
    return None

# Arguments exclusive to SharedLibrary. These are separated to make integrating
# them into build_target easier
_EXCLUSIVE_SHARED_LIB_KWS: T.List[KwargInfo] = [
    _DARWIN_VERSIONS_KW,
    KwargInfo('soversion', (str, int, NoneType), convertor=lambda x: str(x) if x is not None else None),
    KwargInfo('version', (str, NoneType), validator=_validate_shlib_version),
    KwargInfo('shortname', (str, NoneType), since='1.10.0', validator=_shortname_validator),
]

# The total list of arguments used by SharedLibrary
SHARED_LIB_KWS: T.List[KwargInfo] = [
    *_BUILD_TARGET_KWS,
    *_EXCLUSIVE_SHARED_LIB_KWS,
    *_EXCLUSIVE_LIB_KWS,
    _VS_MODULE_DEFS_KW,
    _JAVA_LANG_KW,
]

# Arguments exclusive to SharedModule. These are separated to make integrating
# them into build_target easier
_EXCLUSIVE_SHARED_MOD_KWS: T.List[KwargInfo] = []

# The total list of arguments used by SharedModule
SHARED_MOD_KWS = [
    *_BUILD_TARGET_KWS,
    *_EXCLUSIVE_SHARED_MOD_KWS,
    *_EXCLUSIVE_LIB_KWS,
    _VS_MODULE_DEFS_KW,
    _JAVA_LANG_KW,
]

# Arguments exclusive to JAR. These are separated to make integrating
# them into build_target easier
_EXCLUSIVE_JAR_KWS: T.List[KwargInfo] = [
    KwargInfo('main_class', str, default=''),
    KwargInfo('java_resources', (StructuredSources, NoneType), since='0.62.0'),
    _JAVA_LANG_KW.evolve(deprecated=None, deprecated_message=None),
]

# The total list of arguments used by JAR
JAR_KWS = [
    *_ALL_TARGET_KWS,
    *_EXCLUSIVE_JAR_KWS,
    KwargInfo(
        'sources',
        ContainerTypeInfo(list, (str, File, CustomTarget, CustomTargetIndex, GeneratedList, ExtractedObjects, BuildTarget)),
        listify=True,
        default=[],
    ),
    *[a.evolve(deprecated='1.3.0', deprecated_message='This argument has never done anything in jar(), and should be removed')
      for a in _LANGUAGE_KWS],
]

_SHARED_STATIC_ARGS: T.List[KwargInfo[T.List[str]]] = [
    *[l.evolve(name=l.name.replace('_', '_static_'), since='1.3.0')
      for l in _LANGUAGE_KWS],
    *[l.evolve(name=l.name.replace('_', '_shared_'), since='1.3.0')
      for l in _LANGUAGE_KWS],
]

# Arguments used by both_library and library
LIBRARY_KWS = [
    *_BUILD_TARGET_KWS,
    *_EXCLUSIVE_LIB_KWS,
    *_EXCLUSIVE_SHARED_LIB_KWS,
    *_EXCLUSIVE_SHARED_MOD_KWS,
    *_EXCLUSIVE_STATIC_LIB_KWS,
    *_SHARED_STATIC_ARGS,
    _VS_MODULE_DEFS_KW,
    _JAVA_LANG_KW,
]

# Arguments used by build_Target
BUILD_TARGET_KWS = [
    *_BUILD_TARGET_KWS,
    *_EXCLUSIVE_SHARED_LIB_KWS,
    *_EXCLUSIVE_SHARED_MOD_KWS,
    *_EXCLUSIVE_STATIC_LIB_KWS,
    *EXCLUSIVE_EXECUTABLE_KWS,
    *_SHARED_STATIC_ARGS,
    RUST_ABI_KW.evolve(since='1.10.0'),
    *[a.evolve(deprecated='1.3.0', deprecated_message='The use of "jar" in "build_target()" is deprecated, and this argument is only used by jar()')
      for a in _EXCLUSIVE_JAR_KWS],
    KwargInfo(
        'target_type',
        str,
        required=True,
        validator=in_set_validator({
            'executable', 'shared_library', 'static_library', 'shared_module',
            'both_libraries', 'library', 'jar'
        }),
        since_values={
            'shared_module': '0.51.0',
        },
        deprecated_values={
            'jar': ('1.3.0', 'use the "jar()" function directly'),
        }
    )
]

def _pkgconfig_define_convertor(x: T.List[str]) -> PkgConfigDefineType:
    if x:
        keys = itertools.islice(x, 0, None, 2)
        vals = itertools.islice(x, 1, None, 2)
        return tuple(zip(keys, vals))
    return None

PKGCONFIG_DEFINE_KW: KwargInfo = KwargInfo(
    'pkgconfig_define',
    ContainerTypeInfo(list, str, pairs=True),
    default=[],
    convertor=_pkgconfig_define_convertor,
)

INCLUDE_TYPE = KwargInfo(
    'include_type',
    str,
    default='preserve',
    since='0.52.0',
    validator=in_set_validator({'system', 'non-system', 'preserve'})
)


_DEPRECATED_DEPENDENCY_METHODS = frozenset(
    {'sdlconfig', 'cups-config', 'pcap-config', 'libwmf-config', 'qmake'})


def _dependency_method_convertor(value: str) -> DependencyMethods:
    if value in _DEPRECATED_DEPENDENCY_METHODS:
        return DependencyMethods.CONFIG_TOOL
    return DependencyMethods(value)


DEPENDENCY_METHOD_KW = KwargInfo(
    'method',
    str,
    default='auto',
    since='0.40.0',
    validator=in_set_validator(
        {m.value for m in DependencyMethods} | _DEPRECATED_DEPENDENCY_METHODS),
    convertor=_dependency_method_convertor,
    deprecated_values={
        'sdlconfig': ('0.44.0', 'use config-tool instead'),
        'cups-config': ('0.44.0', 'use config-tool instead'),
        'pcap-config': ('0.44.0', 'use config-tool instead'),
        'libwmf-config': ('0.44.0', 'use config-tool instead'),
        'qmake': ('0.58.0', 'use config-tool instead'),
    },
)


DEPENDENCY_KWS: T.List[KwargInfo] = [
    DEFAULT_OPTIONS.evolve(since='0.38.0'),
    DEPENDENCY_METHOD_KW,
    DISABLER_KW.evolve(since='0.49.0'),
    INCLUDE_TYPE,
    NATIVE_KW,
    REQUIRED_KW,
    KwargInfo('allow_fallback', (bool, NoneType), since='0.56.0'),
    KwargInfo('cmake_args', ContainerTypeInfo(list, str), listify=True, default=[], since='0.50.0'),
    KwargInfo('cmake_module_path', ContainerTypeInfo(list, str), listify=True, default=[], since='0.50.0'),
    KwargInfo('cmake_package_version', str, default='', since='0.57.0'),
    KwargInfo('components', ContainerTypeInfo(list, str), listify=True, default=[], since='0.54.0'),
    KwargInfo('fallback', (ContainerTypeInfo(list, str), str, NoneType), since='0.54.0'),
    KwargInfo('language', (str, NoneType), convertor=lambda x: x.lower() if x is not None else x,
              validator=lambda x: 'Must be a valid language if set' if (x is not None and x not in compilers.all_languages) else None),
    KwargInfo('main', bool, default=False),
    KwargInfo('modules', ContainerTypeInfo(list, str), listify=True, default=[]),
    KwargInfo('not_found_message', str, default='', since='0.50.0'),
    KwargInfo('optional_modules', ContainerTypeInfo(list, str), listify=True, default=[]),
    KwargInfo('private_headers', bool, default=False),
    KwargInfo('static', (bool, NoneType)),
    KwargInfo('version', ContainerTypeInfo(list, str), listify=True, default=[]),
]
