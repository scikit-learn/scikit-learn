from __future__ import annotations
import os
import shlex
import subprocess
import copy
import textwrap

from pathlib import Path, PurePath

from .. import mesonlib
from .. import options
from .. import build
from .. import mlog

from ..modules import ModuleReturnValue, ModuleObject, ModuleState, ExtensionModule, NewExtensionModule
from ..backend.backends import TestProtocol
from ..interpreterbase import (
                               ContainerTypeInfo, KwargInfo, InterpreterObject, MesonOperator,
                               MesonInterpreterObject, ObjectHolder, MutableInterpreterObject,
                               FeatureNew, FeatureDeprecated,
                               typed_pos_args, typed_kwargs, typed_operator,
                               noArgsFlattening, noPosargs, noKwargs, unholder_return,
                               flatten, resolve_second_level_holders, InterpreterException, InvalidArguments, InvalidCode)
from ..interpreter.type_checking import NoneType, ENV_KW, ENV_SEPARATOR_KW, PKGCONFIG_DEFINE_KW
from ..dependencies import Dependency, ExternalLibrary, InternalDependency
from ..programs import ExternalProgram
from ..mesonlib import HoldableObject, listify, Popen_safe

import typing as T

if T.TYPE_CHECKING:
    from . import kwargs
    from ..cmake.interpreter import CMakeInterpreter
    from ..envconfig import MachineInfo
    from ..interpreterbase import FeatureCheckBase, SubProject, TYPE_var, TYPE_kwargs, TYPE_nvar, TYPE_nkwargs
    from .interpreter import Interpreter

    from typing_extensions import TypedDict

    class EnvironmentSeparatorKW(TypedDict):

        separator: str

    class InternalDependencyAsKW(TypedDict):

        recursive: bool

_ERROR_MSG_KW: KwargInfo[T.Optional[str]] = KwargInfo('error_message', (str, NoneType))


def extract_required_kwarg(kwargs: 'kwargs.ExtractRequired',
                           subproject: 'SubProject',
                           feature_check: T.Optional[FeatureCheckBase] = None,
                           default: bool = True) -> T.Tuple[bool, bool, T.Optional[str]]:
    val = kwargs.get('required', default)
    disabled = False
    required = False
    feature: T.Optional[str] = None
    if isinstance(val, options.UserFeatureOption):
        if not feature_check:
            feature_check = FeatureNew('User option "feature"', '0.47.0')
        feature_check.use(subproject)
        feature = val.name
        if val.is_disabled():
            disabled = True
        elif val.is_enabled():
            required = True
    elif isinstance(val, bool):
        required = val
    else:
        raise InterpreterException('required keyword argument must be boolean or a feature option')

    # Keep boolean value in kwargs to simplify other places where this kwarg is
    # checked.
    # TODO: this should be removed, and those callers should learn about FeatureOptions
    kwargs['required'] = required

    return disabled, required, feature

def extract_search_dirs(kwargs: 'kwargs.ExtractSearchDirs') -> T.List[str]:
    search_dirs_str = mesonlib.stringlistify(kwargs.get('dirs', []))
    search_dirs = [Path(d).expanduser() for d in search_dirs_str]
    for d in search_dirs:
        if mesonlib.is_windows() and d.root.startswith('\\'):
            # a Unix-path starting with `/` that is not absolute on Windows.
            # discard without failing for end-user ease of cross-platform directory arrays
            continue
        if not d.is_absolute():
            raise InvalidCode(f'Search directory {d} is not an absolute path.')
    return [str(s) for s in search_dirs]

class FeatureOptionHolder(ObjectHolder[options.UserFeatureOption]):
    def __init__(self, option: options.UserFeatureOption, interpreter: 'Interpreter'):
        super().__init__(option, interpreter)
        if option and option.is_auto():
            # TODO: we need to cast here because options is not a TypedDict
            auto = T.cast('options.UserFeatureOption', self.env.coredata.optstore.get_value_object_for('auto_features'))
            self.held_object = copy.copy(auto)
            self.held_object.name = option.name

    @property
    def value(self) -> str:
        return 'disabled' if not self.held_object else self.held_object.value

    def as_disabled(self) -> options.UserFeatureOption:
        disabled = copy.deepcopy(self.held_object)
        disabled.value = 'disabled'
        return disabled

    def as_enabled(self) -> options.UserFeatureOption:
        enabled = copy.deepcopy(self.held_object)
        enabled.value = 'enabled'
        return enabled

    @noPosargs
    @noKwargs
    @InterpreterObject.method('enabled')
    def enabled_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> bool:
        return self.value == 'enabled'

    @noPosargs
    @noKwargs
    @InterpreterObject.method('disabled')
    def disabled_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> bool:
        return self.value == 'disabled'

    @noPosargs
    @noKwargs
    @FeatureNew('feature_option.allowed()', '0.59.0')
    @InterpreterObject.method('allowed')
    def allowed_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> bool:
        return self.value != 'disabled'

    @noPosargs
    @noKwargs
    @InterpreterObject.method('auto')
    def auto_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> bool:
        return self.value == 'auto'

    def _disable_if(self, condition: bool, message: T.Optional[str]) -> options.UserFeatureOption:
        if not condition:
            return copy.deepcopy(self.held_object)

        if self.value == 'enabled':
            err_msg = f'Feature {self.held_object.name} cannot be enabled'
            if message:
                err_msg += f': {message}'
            raise InterpreterException(err_msg)
        return self.as_disabled()

    @FeatureNew('feature_option.require()', '0.59.0')
    @typed_pos_args('feature_option.require', bool)
    @typed_kwargs(
        'feature_option.require',
        _ERROR_MSG_KW,
    )
    @InterpreterObject.method('require')
    def require_method(self, args: T.Tuple[bool], kwargs: 'kwargs.FeatureOptionRequire') -> options.UserFeatureOption:
        return self._disable_if(not args[0], kwargs['error_message'])

    @FeatureNew('feature_option.disable_if()', '1.1.0')
    @typed_pos_args('feature_option.disable_if', bool)
    @typed_kwargs(
        'feature_option.disable_if',
        _ERROR_MSG_KW,
    )
    @InterpreterObject.method('disable_if')
    def disable_if_method(self, args: T.Tuple[bool], kwargs: 'kwargs.FeatureOptionRequire') -> options.UserFeatureOption:
        return self._disable_if(args[0], kwargs['error_message'])

    @FeatureNew('feature_option.enable_if()', '1.1.0')
    @typed_pos_args('feature_option.enable_if', bool)
    @typed_kwargs(
        'feature_option.enable_if',
        _ERROR_MSG_KW,
    )
    @InterpreterObject.method('enable_if')
    def enable_if_method(self, args: T.Tuple[bool], kwargs: 'kwargs.FeatureOptionRequire') -> options.UserFeatureOption:
        if not args[0]:
            return copy.deepcopy(self.held_object)

        if self.value == 'disabled':
            err_msg = f'Feature {self.held_object.name} cannot be disabled'
            if kwargs['error_message']:
                err_msg += f': {kwargs["error_message"]}'
            raise InterpreterException(err_msg)
        return self.as_enabled()

    @FeatureNew('feature_option.disable_auto_if()', '0.59.0')
    @noKwargs
    @typed_pos_args('feature_option.disable_auto_if', bool)
    @InterpreterObject.method('disable_auto_if')
    def disable_auto_if_method(self, args: T.Tuple[bool], kwargs: TYPE_kwargs) -> options.UserFeatureOption:
        return copy.deepcopy(self.held_object) if self.value != 'auto' or not args[0] else self.as_disabled()

    @FeatureNew('feature_option.enable_auto_if()', '1.1.0')
    @noKwargs
    @typed_pos_args('feature_option.enable_auto_if', bool)
    @InterpreterObject.method('enable_auto_if')
    def enable_auto_if_method(self, args: T.Tuple[bool], kwargs: TYPE_kwargs) -> options.UserFeatureOption:
        return self.as_enabled() if self.value == 'auto' and args[0] else copy.deepcopy(self.held_object)


class RunProcess(MesonInterpreterObject):

    def __init__(self,
                 cmd: ExternalProgram,
                 args: T.List[str],
                 env: mesonlib.EnvironmentVariables,
                 source_dir: str,
                 build_dir: str,
                 subdir: str,
                 mesonintrospect: T.List[str],
                 in_builddir: bool = False,
                 check: bool = False,
                 capture: bool = True) -> None:
        super().__init__()
        if not isinstance(cmd, ExternalProgram):
            raise AssertionError('BUG: RunProcess must be passed an ExternalProgram')
        self.capture = capture
        self.returncode, self.stdout, self.stderr = self.run_command(cmd, args, env, source_dir, build_dir, subdir, mesonintrospect, in_builddir, check)

    def run_command(self,
                    cmd: ExternalProgram,
                    args: T.List[str],
                    env: mesonlib.EnvironmentVariables,
                    source_dir: str,
                    build_dir: str,
                    subdir: str,
                    mesonintrospect: T.List[str],
                    in_builddir: bool,
                    check: bool = False) -> T.Tuple[int, str, str]:
        command_array = cmd.get_command() + args
        menv = {'MESON_SOURCE_ROOT': source_dir,
                'MESON_BUILD_ROOT': build_dir,
                'MESON_SUBDIR': subdir,
                'MESONINTROSPECT': ' '.join([shlex.quote(x) for x in mesonintrospect]),
                }
        if in_builddir:
            cwd = os.path.join(build_dir, subdir)
        else:
            cwd = os.path.join(source_dir, subdir)
        child_env = os.environ.copy()
        child_env.update(menv)
        child_env = env.get_env(child_env)
        stdout = subprocess.PIPE if self.capture else subprocess.DEVNULL
        mlog.debug('Running command:', mesonlib.join_args(command_array))
        try:
            p, o, e = Popen_safe(command_array, stdout=stdout, env=child_env, cwd=cwd)
            if self.capture:
                mlog.debug('--- stdout ---')
                mlog.debug(o)
            else:
                o = ''
                mlog.debug('--- stdout disabled ---')
            mlog.debug('--- stderr ---')
            mlog.debug(e)
            mlog.debug('')

            if check and p.returncode != 0:
                raise InterpreterException('Command `{}` failed with status {}.'.format(mesonlib.join_args(command_array), p.returncode))

            return p.returncode, o, e
        except FileNotFoundError:
            raise InterpreterException('Could not execute command `%s`.' % mesonlib.join_args(command_array))

    @noPosargs
    @noKwargs
    @InterpreterObject.method('returncode')
    def returncode_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> int:
        return self.returncode

    @noPosargs
    @noKwargs
    @InterpreterObject.method('stdout')
    def stdout_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self.stdout

    @noPosargs
    @noKwargs
    @InterpreterObject.method('stderr')
    def stderr_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self.stderr

class EnvironmentVariablesHolder(ObjectHolder[mesonlib.EnvironmentVariables], MutableInterpreterObject):

    def __init__(self, obj: mesonlib.EnvironmentVariables, interpreter: 'Interpreter'):
        super().__init__(obj, interpreter)

    def __repr__(self) -> str:
        repr_str = "<{0}: {1}>"
        return repr_str.format(self.__class__.__name__, self.held_object.envvars)

    def __deepcopy__(self, memo: T.Dict[str, object]) -> 'EnvironmentVariablesHolder':
        # Avoid trying to copy the interpreter
        return EnvironmentVariablesHolder(copy.deepcopy(self.held_object), self.interpreter)

    def warn_if_has_name(self, name: str) -> None:
        # Multiple append/prepend operations was not supported until 0.58.0.
        if self.held_object.has_name(name):
            m = f'Overriding previous value of environment variable {name!r} with a new one'
            FeatureNew(m, '0.58.0').use(self.subproject, self.current_node)

    @typed_pos_args('environment.set', str, varargs=str, min_varargs=1)
    @typed_kwargs('environment.set', ENV_SEPARATOR_KW)
    @InterpreterObject.method('set')
    def set_method(self, args: T.Tuple[str, T.List[str]], kwargs: 'EnvironmentSeparatorKW') -> None:
        name, values = args
        self.held_object.set(name, values, kwargs['separator'])

    @FeatureNew('environment.unset', '1.4.0')
    @typed_pos_args('environment.unset', str)
    @noKwargs
    @InterpreterObject.method('unset')
    def unset_method(self, args: T.Tuple[str], kwargs: TYPE_kwargs) -> None:
        self.held_object.unset(args[0])

    @typed_pos_args('environment.append', str, varargs=str, min_varargs=1)
    @typed_kwargs('environment.append', ENV_SEPARATOR_KW)
    @InterpreterObject.method('append')
    def append_method(self, args: T.Tuple[str, T.List[str]], kwargs: 'EnvironmentSeparatorKW') -> None:
        name, values = args
        self.warn_if_has_name(name)
        self.held_object.append(name, values, kwargs['separator'])

    @typed_pos_args('environment.prepend', str, varargs=str, min_varargs=1)
    @typed_kwargs('environment.prepend', ENV_SEPARATOR_KW)
    @InterpreterObject.method('prepend')
    def prepend_method(self, args: T.Tuple[str, T.List[str]], kwargs: 'EnvironmentSeparatorKW') -> None:
        name, values = args
        self.warn_if_has_name(name)
        self.held_object.prepend(name, values, kwargs['separator'])


_CONF_DATA_SET_KWS: KwargInfo[T.Optional[str]] = KwargInfo('description', (str, NoneType))


class ConfigurationDataHolder(ObjectHolder[build.ConfigurationData], MutableInterpreterObject):

    def __init__(self, obj: build.ConfigurationData, interpreter: 'Interpreter'):
        super().__init__(obj, interpreter)

    def __deepcopy__(self, memo: T.Dict) -> 'ConfigurationDataHolder':
        return ConfigurationDataHolder(copy.deepcopy(self.held_object), self.interpreter)

    def is_used(self) -> bool:
        return self.held_object.used

    def __check_used(self) -> None:
        if self.is_used():
            raise InterpreterException("Can not set values on configuration object that has been used.")

    @typed_pos_args('configuration_data.set', str, (str, int, bool))
    @typed_kwargs('configuration_data.set', _CONF_DATA_SET_KWS)
    @InterpreterObject.method('set')
    def set_method(self, args: T.Tuple[str, T.Union[str, int, bool]], kwargs: 'kwargs.ConfigurationDataSet') -> None:
        self.__check_used()
        self.held_object.values[args[0]] = (args[1], kwargs['description'])

    @typed_pos_args('configuration_data.set_quoted', str, str)
    @typed_kwargs('configuration_data.set_quoted', _CONF_DATA_SET_KWS)
    @InterpreterObject.method('set_quoted')
    def set_quoted_method(self, args: T.Tuple[str, str], kwargs: 'kwargs.ConfigurationDataSet') -> None:
        self.__check_used()
        escaped_val = '\\"'.join(args[1].split('"'))
        self.held_object.values[args[0]] = (f'"{escaped_val}"', kwargs['description'])

    @typed_pos_args('configuration_data.set10', str, (int, bool))
    @typed_kwargs('configuration_data.set10', _CONF_DATA_SET_KWS)
    @InterpreterObject.method('set10')
    def set10_method(self, args: T.Tuple[str, T.Union[int, bool]], kwargs: 'kwargs.ConfigurationDataSet') -> None:
        self.__check_used()
        # bool is a subclass of int, so we need to check for bool explicitly.
        # We already have typed_pos_args checking that this is either a bool or
        # an int.
        if not isinstance(args[1], bool):
            mlog.deprecation('configuration_data.set10 with number. The `set10` '
                             'method should only be used with booleans',
                             location=self.interpreter.current_node)
            if args[1] < 0:
                mlog.warning('Passing a number that is less than 0 may not have the intended result, '
                             'as meson will treat all non-zero values as true.',
                             location=self.interpreter.current_node)
        self.held_object.values[args[0]] = (int(args[1]), kwargs['description'])

    @typed_pos_args('configuration_data.has', (str, int, bool))
    @noKwargs
    @InterpreterObject.method('has')
    def has_method(self, args: T.Tuple[T.Union[str, int, bool]], kwargs: TYPE_kwargs) -> bool:
        return args[0] in self.held_object.values

    @FeatureNew('configuration_data.get()', '0.38.0')
    @typed_pos_args('configuration_data.get', str, optargs=[(str, int, bool)])
    @noKwargs
    @InterpreterObject.method('get')
    def get_method(self, args: T.Tuple[str, T.Optional[T.Union[str, int, bool]]],
                   kwargs: TYPE_kwargs) -> T.Union[str, int, bool]:
        name = args[0]
        if name in self.held_object:
            return self.held_object.get(name)[0]
        elif args[1] is not None:
            return args[1]
        raise InterpreterException(f'Entry {name} not in configuration data.')

    @FeatureNew('configuration_data.get_unquoted()', '0.44.0')
    @typed_pos_args('configuration_data.get_unquoted', str, optargs=[(str, int, bool)])
    @noKwargs
    @InterpreterObject.method('get_unquoted')
    def get_unquoted_method(self, args: T.Tuple[str, T.Optional[T.Union[str, int, bool]]],
                            kwargs: TYPE_kwargs) -> T.Union[str, int, bool]:
        name = args[0]
        if name in self.held_object:
            val = self.held_object.get(name)[0]
        elif args[1] is not None:
            val = args[1]
        else:
            raise InterpreterException(f'Entry {name} not in configuration data.')
        if isinstance(val, str) and val[0] == '"' and val[-1] == '"':
            return val[1:-1]
        return val

    def get(self, name: str) -> T.Tuple[T.Union[str, int, bool], T.Optional[str]]:
        return self.held_object.values[name]

    @FeatureNew('configuration_data.keys()', '0.57.0')
    @noPosargs
    @noKwargs
    @InterpreterObject.method('keys')
    def keys_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> T.List[str]:
        return sorted(self.keys())

    def keys(self) -> T.List[str]:
        return list(self.held_object.values.keys())

    @typed_pos_args('configuration_data.merge_from', build.ConfigurationData)
    @noKwargs
    @InterpreterObject.method('merge_from')
    def merge_from_method(self, args: T.Tuple[build.ConfigurationData], kwargs: TYPE_kwargs) -> None:
        from_object = args[0]
        self.held_object.values.update(from_object.values)


_PARTIAL_DEP_KWARGS = [
    KwargInfo('compile_args', bool, default=False),
    KwargInfo('link_args',    bool, default=False),
    KwargInfo('links',        bool, default=False),
    KwargInfo('includes',     bool, default=False),
    KwargInfo('sources',      bool, default=False),
]

class DependencyHolder(ObjectHolder[Dependency]):
    def __init__(self, dep: Dependency, interpreter: 'Interpreter'):
        super().__init__(dep, interpreter)

    def found(self) -> bool:
        return self.found_method([], {})

    @noPosargs
    @noKwargs
    @InterpreterObject.method('type_name')
    def type_name_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self.held_object.type_name

    @noPosargs
    @noKwargs
    @InterpreterObject.method('found')
    def found_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> bool:
        if self.held_object.type_name == 'internal':
            return True
        return self.held_object.found()

    @noPosargs
    @noKwargs
    @InterpreterObject.method('version')
    def version_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self.held_object.get_version()

    @noPosargs
    @noKwargs
    @InterpreterObject.method('name')
    def name_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self.held_object.get_name()

    @FeatureDeprecated('dependency.get_pkgconfig_variable', '0.56.0',
                       'use dependency.get_variable(pkgconfig : ...) instead')
    @typed_pos_args('dependency.get_pkgconfig_variable', str)
    @typed_kwargs(
        'dependency.get_pkgconfig_variable',
        KwargInfo('default', str, default=''),
        PKGCONFIG_DEFINE_KW.evolve(name='define_variable')
    )
    @InterpreterObject.method('get_pkgconfig_variable')
    def pkgconfig_method(self, args: T.Tuple[str], kwargs: 'kwargs.DependencyPkgConfigVar') -> str:
        from ..dependencies.pkgconfig import PkgConfigDependency
        if not isinstance(self.held_object, PkgConfigDependency):
            raise InvalidArguments(f'{self.held_object.get_name()!r} is not a pkgconfig dependency')
        if kwargs['define_variable'] and len(kwargs['define_variable']) > 1:
            FeatureNew.single_use('dependency.get_pkgconfig_variable keyword argument "define_variable"  with more than one pair',
                                  '1.3.0', self.subproject, location=self.current_node)
        return self.held_object.get_variable(
            pkgconfig=args[0],
            default_value=kwargs['default'],
            pkgconfig_define=kwargs['define_variable'],
        )

    @FeatureNew('dependency.get_configtool_variable', '0.44.0')
    @FeatureDeprecated('dependency.get_configtool_variable', '0.56.0',
                       'use dependency.get_variable(configtool : ...) instead')
    @noKwargs
    @typed_pos_args('dependency.get_config_tool_variable', str)
    @InterpreterObject.method('get_configtool_variable')
    def configtool_method(self, args: T.Tuple[str], kwargs: TYPE_kwargs) -> str:
        from ..dependencies.configtool import ConfigToolDependency
        if not isinstance(self.held_object, ConfigToolDependency):
            raise InvalidArguments(f'{self.held_object.get_name()!r} is not a config-tool dependency')
        return self.held_object.get_variable(
            configtool=args[0],
            default_value='',
        )

    @FeatureNew('dependency.partial_dependency', '0.46.0')
    @noPosargs
    @typed_kwargs('dependency.partial_dependency', *_PARTIAL_DEP_KWARGS)
    @InterpreterObject.method('partial_dependency')
    def partial_dependency_method(self, args: T.List[TYPE_nvar], kwargs: 'kwargs.DependencyMethodPartialDependency') -> Dependency:
        pdep = self.held_object.get_partial_dependency(**kwargs)
        return pdep

    @FeatureNew('dependency.get_variable', '0.51.0')
    @typed_pos_args('dependency.get_variable', optargs=[str])
    @typed_kwargs(
        'dependency.get_variable',
        KwargInfo('cmake', (str, NoneType)),
        KwargInfo('pkgconfig', (str, NoneType)),
        KwargInfo('configtool', (str, NoneType)),
        KwargInfo('internal', (str, NoneType), since='0.54.0'),
        KwargInfo('system', (str, NoneType), since='1.6.0'),
        KwargInfo('default_value', (str, NoneType)),
        PKGCONFIG_DEFINE_KW,
    )
    @InterpreterObject.method('get_variable')
    def variable_method(self, args: T.Tuple[T.Optional[str]], kwargs: 'kwargs.DependencyGetVariable') -> str:
        default_varname = args[0]
        if default_varname is not None:
            FeatureNew('Positional argument to dependency.get_variable()', '0.58.0').use(self.subproject, self.current_node)
        if kwargs['pkgconfig_define'] and len(kwargs['pkgconfig_define']) > 1:
            FeatureNew.single_use('dependency.get_variable keyword argument "pkgconfig_define" with more than one pair',
                                  '1.3.0', self.subproject, 'In previous versions, this silently returned a malformed value.',
                                  self.current_node)
        return self.held_object.get_variable(
            cmake=kwargs['cmake'] or default_varname,
            pkgconfig=kwargs['pkgconfig'] or default_varname,
            configtool=kwargs['configtool'] or default_varname,
            internal=kwargs['internal'] or default_varname,
            system=kwargs['system'] or default_varname,
            default_value=kwargs['default_value'],
            pkgconfig_define=kwargs['pkgconfig_define'],
        )

    @FeatureNew('dependency.include_type', '0.52.0')
    @noPosargs
    @noKwargs
    @InterpreterObject.method('include_type')
    def include_type_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self.held_object.get_include_type()

    @FeatureNew('dependency.as_system', '0.52.0')
    @noKwargs
    @typed_pos_args('dependency.as_system', optargs=[str])
    @InterpreterObject.method('as_system')
    def as_system_method(self, args: T.Tuple[T.Optional[str]], kwargs: TYPE_kwargs) -> Dependency:
        return self.held_object.generate_system_dependency(args[0] or 'system')

    @FeatureNew('dependency.as_link_whole', '0.56.0')
    @noKwargs
    @noPosargs
    @InterpreterObject.method('as_link_whole')
    def as_link_whole_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> Dependency:
        if not isinstance(self.held_object, InternalDependency):
            raise InterpreterException('as_link_whole method is only supported on declare_dependency() objects')
        new_dep = self.held_object.generate_link_whole_dependency()
        return new_dep

    @FeatureNew('dependency.as_static', '1.6.0')
    @noPosargs
    @typed_kwargs(
        'dependency.as_static',
        KwargInfo('recursive', bool, default=False),
    )
    @InterpreterObject.method('as_static')
    def as_static_method(self, args: T.List[TYPE_var], kwargs: InternalDependencyAsKW) -> Dependency:
        if not isinstance(self.held_object, InternalDependency):
            raise InterpreterException('as_static method is only supported on declare_dependency() objects')
        return self.held_object.get_as_static(kwargs['recursive'])

    @FeatureNew('dependency.as_shared', '1.6.0')
    @noPosargs
    @typed_kwargs(
        'dependency.as_shared',
        KwargInfo('recursive', bool, default=False),
    )
    @InterpreterObject.method('as_shared')
    def as_shared_method(self, args: T.List[TYPE_var], kwargs: InternalDependencyAsKW) -> Dependency:
        if not isinstance(self.held_object, InternalDependency):
            raise InterpreterException('as_shared method is only supported on declare_dependency() objects')
        return self.held_object.get_as_shared(kwargs['recursive'])

_EXTPROG = T.TypeVar('_EXTPROG', bound=ExternalProgram)

class _ExternalProgramHolder(ObjectHolder[_EXTPROG]):
    def __init__(self, ep: _EXTPROG, interpreter: 'Interpreter') -> None:
        super().__init__(ep, interpreter)

    @noPosargs
    @noKwargs
    @InterpreterObject.method('found')
    def found_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> bool:
        return self.found()

    @noPosargs
    @noKwargs
    @FeatureDeprecated('ExternalProgram.path', '0.55.0',
                       'use ExternalProgram.full_path() instead')
    @InterpreterObject.method('path')
    def path_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self._full_path()

    @noPosargs
    @noKwargs
    @FeatureNew('ExternalProgram.full_path', '0.55.0')
    @InterpreterObject.method('full_path')
    def full_path_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self._full_path()

    def _full_path(self) -> str:
        if not self.found():
            raise InterpreterException('Unable to get the path of a not-found external program')
        path = self.held_object.get_path()
        assert path is not None
        return path

    @noPosargs
    @noKwargs
    @FeatureNew('ExternalProgram.version', '0.62.0')
    @InterpreterObject.method('version')
    def version_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        if not self.found():
            raise InterpreterException('Unable to get the version of a not-found external program')
        try:
            return self.held_object.get_version(self.interpreter)
        except mesonlib.MesonException:
            return 'unknown'

    def found(self) -> bool:
        return self.held_object.found()

class ExternalProgramHolder(_ExternalProgramHolder[ExternalProgram]):
    pass

class ExternalLibraryHolder(ObjectHolder[ExternalLibrary]):
    def __init__(self, el: ExternalLibrary, interpreter: 'Interpreter'):
        super().__init__(el, interpreter)

    @noPosargs
    @noKwargs
    @InterpreterObject.method('type_name')
    def type_name_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self.held_object.type_name

    @noPosargs
    @noKwargs
    @InterpreterObject.method('found')
    def found_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> bool:
        return self.held_object.found()

    @FeatureNew('dependency.partial_dependency', '0.46.0')
    @noPosargs
    @typed_kwargs('dependency.partial_dependency', *_PARTIAL_DEP_KWARGS)
    @InterpreterObject.method('partial_dependency')
    def partial_dependency_method(self, args: T.List[TYPE_nvar], kwargs: 'kwargs.DependencyMethodPartialDependency') -> Dependency:
        pdep = self.held_object.get_partial_dependency(**kwargs)
        return pdep

    @FeatureNew('dependency.name', '1.5.0')
    @noPosargs
    @noKwargs
    @InterpreterObject.method('name')
    def name_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self.held_object.name


# A machine that's statically known from the cross file
class MachineHolder(ObjectHolder['MachineInfo']):
    def __init__(self, machine_info: 'MachineInfo', interpreter: 'Interpreter'):
        super().__init__(machine_info, interpreter)

    @noPosargs
    @noKwargs
    @InterpreterObject.method('cpu_family')
    def cpu_family_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self.held_object.cpu_family

    @noPosargs
    @noKwargs
    @InterpreterObject.method('cpu')
    def cpu_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self.held_object.cpu

    @noPosargs
    @noKwargs
    @InterpreterObject.method('system')
    def system_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self.held_object.system

    @noPosargs
    @noKwargs
    @InterpreterObject.method('endian')
    def endian_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self.held_object.endian

    @noPosargs
    @noKwargs
    @InterpreterObject.method('kernel')
    def kernel_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        if self.held_object.kernel is not None:
            return self.held_object.kernel
        raise InterpreterException('Kernel not defined or could not be autodetected.')

    @noPosargs
    @noKwargs
    @InterpreterObject.method('subsystem')
    def subsystem_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        if self.held_object.subsystem is not None:
            return self.held_object.subsystem
        raise InterpreterException('Subsystem not defined or could not be autodetected.')


class IncludeDirsHolder(ObjectHolder[build.IncludeDirs]):
    pass

class FileHolder(ObjectHolder[mesonlib.File]):
    def __init__(self, file: mesonlib.File, interpreter: 'Interpreter'):
        super().__init__(file, interpreter)

    @noPosargs
    @noKwargs
    @FeatureNew('file.full_path', '1.4.0')
    @InterpreterObject.method('full_path')
    def full_path_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self.held_object.absolute_path(self.env.source_dir, self.env.build_dir)

class HeadersHolder(ObjectHolder[build.Headers]):
    pass

class DataHolder(ObjectHolder[build.Data]):
    pass

class SymlinkDataHolder(ObjectHolder[build.SymlinkData]):
    pass

class InstallDirHolder(ObjectHolder[build.InstallDir]):
    pass

class ManHolder(ObjectHolder[build.Man]):
    pass

class EmptyDirHolder(ObjectHolder[build.EmptyDir]):
    pass

class GeneratedObjectsHolder(ObjectHolder[build.ExtractedObjects]):
    pass

class Test(MesonInterpreterObject):
    def __init__(self, name: str, project: str, suite: T.List[str],
                 exe: T.Union[ExternalProgram, build.Executable, build.CustomTarget, build.CustomTargetIndex],
                 depends: T.List[T.Union[build.CustomTarget, build.BuildTarget]],
                 is_parallel: bool,
                 cmd_args: T.List[T.Union[str, mesonlib.File, build.Target, ExternalProgram]],
                 env: mesonlib.EnvironmentVariables,
                 should_fail: bool, timeout: int, workdir: T.Optional[str], protocol: str,
                 priority: int, verbose: bool):
        super().__init__()
        self.name = name
        self.suite = listify(suite)
        self.project_name = project
        self.exe = exe
        self.depends = depends
        self.is_parallel = is_parallel
        self.cmd_args = cmd_args
        self.env = env
        self.should_fail = should_fail
        self.timeout = timeout
        self.workdir = workdir
        self.protocol = TestProtocol.from_str(protocol)
        self.priority = priority
        self.verbose = verbose

    def get_exe(self) -> T.Union[ExternalProgram, build.Executable, build.CustomTarget, build.CustomTargetIndex]:
        return self.exe

    def get_name(self) -> str:
        return self.name


class Doctest(Test):
    target: T.Optional[build.BuildTarget] = None


class NullSubprojectInterpreter(HoldableObject):
    pass

# TODO: This should really be an `ObjectHolder`, but the additional stuff in this
#       class prevents this. Thus, this class should be split into a pure
#       `ObjectHolder` and a class specifically for storing in `Interpreter`.
class SubprojectHolder(MesonInterpreterObject):

    def __init__(self, subinterpreter: T.Union['Interpreter', NullSubprojectInterpreter],
                 subdir: str,
                 warnings: int = 0,
                 disabled_feature: T.Optional[str] = None,
                 exception: T.Optional[Exception] = None,
                 callstack: T.Optional[T.List[str]] = None) -> None:
        super().__init__()
        self.held_object = subinterpreter
        self.warnings = warnings
        self.disabled_feature = disabled_feature
        self.exception = exception
        self.subdir = PurePath(subdir).as_posix()
        self.cm_interpreter: T.Optional[CMakeInterpreter] = None
        self.callstack = callstack

    @noPosargs
    @noKwargs
    @InterpreterObject.method('found')
    def found_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> bool:
        return self.found()

    def found(self) -> bool:
        return not isinstance(self.held_object, NullSubprojectInterpreter)

    @unholder_return
    def get_variable(self, args: T.Tuple[str, T.Optional[str]], kwargs: TYPE_kwargs) -> T.Union[TYPE_var, InterpreterObject]:
        if isinstance(self.held_object, NullSubprojectInterpreter):  # == not self.found()
            raise InterpreterException(f'Subproject "{self.subdir}" disabled can\'t get_variable on it.')
        varname, fallback = args
        try:
            return self.held_object.variables[varname]
        except KeyError:
            if fallback is not None:
                return self.held_object._holderify(fallback)
            raise InvalidArguments(f'Requested variable "{varname}" not found.')

    @noKwargs
    @typed_pos_args('subproject.get_variable', str, optargs=[object])
    @noArgsFlattening
    @InterpreterObject.method('get_variable')
    def get_variable_method(self, args: T.Tuple[str, T.Optional[str]], kwargs: TYPE_kwargs) -> T.Union[TYPE_var, InterpreterObject]:
        return self.get_variable(args, kwargs)

class ModuleObjectHolder(ObjectHolder[ModuleObject]):
    def method_call(self, method_name: str, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> TYPE_var:
        modobj = self.held_object
        method = modobj.methods.get(method_name)
        if not method:
            raise InvalidCode(f'Unknown method {method_name!r} in object.')
        if not getattr(method, 'no-args-flattening', False):
            args = flatten(args)
        if not getattr(method, 'no-second-level-holder-flattening', False):
            args, kwargs = resolve_second_level_holders(args, kwargs)
        if not self.interpreter.active_projectname:
            assert isinstance(modobj, (ExtensionModule, NewExtensionModule)), 'for mypy'
            full_method_name = f'{modobj.INFO.name}.{method_name}'
            raise mesonlib.MesonException(f'Module methods ({full_method_name}) cannot be invoked during project declaration.')
        state = ModuleState(self.interpreter)
        # Many modules do for example self.interpreter.find_program_impl(),
        # so we have to ensure they use the current interpreter and not the one
        # that first imported that module, otherwise it will use outdated
        # overrides.
        if isinstance(modobj, ExtensionModule):
            modobj.interpreter = self.interpreter
        ret = method(state, args, kwargs)
        if isinstance(ret, ModuleReturnValue):
            self.interpreter.process_new_values(ret.new_objects)
            ret = ret.return_value
        return ret

class MutableModuleObjectHolder(ModuleObjectHolder, MutableInterpreterObject):
    def __deepcopy__(self, memo: T.Dict[int, T.Any]) -> 'MutableModuleObjectHolder':
        # Deepcopy only held object, not interpreter
        modobj = copy.deepcopy(self.held_object, memo)
        return MutableModuleObjectHolder(modobj, self.interpreter)


_BuildTarget = T.TypeVar('_BuildTarget', bound=T.Union[build.BuildTarget, build.BothLibraries])

class BuildTargetHolder(ObjectHolder[_BuildTarget]):
    def __init__(self, target: _BuildTarget, interp: 'Interpreter'):
        super().__init__(target, interp)

    def __repr__(self) -> str:
        r = '<{} {}: {}>'
        h = self.held_object
        assert isinstance(h, build.BuildTarget)
        return r.format(self.__class__.__name__, h.get_id(), h.filename)

    @property
    def _target_object(self) -> build.BuildTarget:
        if isinstance(self.held_object, build.BothLibraries):
            return self.held_object.get_default_object()
        assert isinstance(self.held_object, build.BuildTarget)
        return self.held_object

    def is_cross(self) -> bool:
        return not self._target_object.environment.machines.matches_build_machine(self._target_object.for_machine)

    @noPosargs
    @noKwargs
    @InterpreterObject.method('found')
    def found_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> bool:
        if not (isinstance(self.held_object, build.Executable) and self.held_object.was_returned_by_find_program):
            FeatureNew.single_use('BuildTarget.found', '0.59.0', subproject=self.held_object.subproject)
        return True

    @noPosargs
    @noKwargs
    @InterpreterObject.method('private_dir_include')
    def private_dir_include_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> build.IncludeDirs:
        return build.IncludeDirs('', [], False, [self.interpreter.backend.get_target_private_dir(self._target_object)])

    @noPosargs
    @noKwargs
    @InterpreterObject.method('full_path')
    def full_path_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self.interpreter.backend.get_target_filename_abs(self._target_object)

    @noPosargs
    @noKwargs
    @FeatureDeprecated('BuildTarget.path', '0.55.0', 'Use BuildTarget.full_path instead')
    @InterpreterObject.method('path')
    def path_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self.interpreter.backend.get_target_filename_abs(self._target_object)

    @noPosargs
    @noKwargs
    @InterpreterObject.method('outdir')
    def outdir_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self.interpreter.backend.get_target_dir(self._target_object)

    @noKwargs
    @typed_pos_args('extract_objects', varargs=(mesonlib.File, str, build.CustomTarget, build.CustomTargetIndex, build.GeneratedList))
    @InterpreterObject.method('extract_objects')
    def extract_objects_method(self, args: T.Tuple[T.List[T.Union[mesonlib.FileOrString, 'build.GeneratedTypes']]], kwargs: TYPE_nkwargs) -> build.ExtractedObjects:
        tobj = self._target_object
        unity_value = self.interpreter.coredata.get_option_for_target(tobj, "unity")
        is_unity = (unity_value == 'on' or (unity_value == 'subprojects' and tobj.subproject != ''))
        return tobj.extract_objects(args[0], is_unity)

    @noPosargs
    @typed_kwargs(
        'extract_all_objects',
        KwargInfo(
            'recursive', bool, default=False, since='0.46.0',
            not_set_warning=textwrap.dedent('''\
                extract_all_objects called without setting recursive
                keyword argument. Meson currently defaults to
                non-recursive to maintain backward compatibility but
                the default will be changed in meson 2.0.
            ''')
        )
    )
    @InterpreterObject.method('extract_all_objects')
    def extract_all_objects_method(self, args: T.List[TYPE_nvar], kwargs: 'kwargs.BuildTargeMethodExtractAllObjects') -> build.ExtractedObjects:
        return self._target_object.extract_all_objects(kwargs['recursive'])

    @noPosargs
    @noKwargs
    @FeatureDeprecated('BuildTarget.get_id', '1.2.0',
                       'This was never formally documented and does not seem to have a real world use. ' +
                       'See https://github.com/mesonbuild/meson/pull/6061')
    @InterpreterObject.method('get_id')
    def get_id_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self._target_object.get_id()

    @FeatureNew('name', '0.54.0')
    @noPosargs
    @noKwargs
    @InterpreterObject.method('name')
    def name_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self._target_object.name

class ExecutableHolder(BuildTargetHolder[build.Executable]):
    pass

class StaticLibraryHolder(BuildTargetHolder[build.StaticLibrary]):
    pass

class SharedLibraryHolder(BuildTargetHolder[build.SharedLibrary]):
    pass

class BothLibrariesHolder(BuildTargetHolder[build.BothLibraries]):
    def __init__(self, libs: build.BothLibraries, interp: 'Interpreter'):
        super().__init__(libs, interp)

    def __repr__(self) -> str:
        r = '<{} {}: {}, {}: {}>'
        h1 = self.held_object.shared
        h2 = self.held_object.static
        return r.format(self.__class__.__name__, h1.get_id(), h1.filename, h2.get_id(), h2.filename)

    @noPosargs
    @noKwargs
    @InterpreterObject.method('get_shared_lib')
    def get_shared_lib_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> build.SharedLibrary:
        lib = copy.copy(self.held_object.shared)
        lib.both_lib = None
        return lib

    @noPosargs
    @noKwargs
    @InterpreterObject.method('get_static_lib')
    def get_static_lib_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> build.StaticLibrary:
        lib = copy.copy(self.held_object.static)
        lib.both_lib = None
        return lib

class SharedModuleHolder(BuildTargetHolder[build.SharedModule]):
    pass

class JarHolder(BuildTargetHolder[build.Jar]):
    pass

class CustomTargetIndexHolder(ObjectHolder[build.CustomTargetIndex]):
    def __init__(self, target: build.CustomTargetIndex, interp: 'Interpreter'):
        super().__init__(target, interp)

    @FeatureNew('custom_target[i].full_path', '0.54.0')
    @noPosargs
    @noKwargs
    @InterpreterObject.method('full_path')
    def full_path_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        assert self.interpreter.backend is not None
        return self.interpreter.backend.get_target_filename_abs(self.held_object)

_CT = T.TypeVar('_CT', bound=build.CustomTarget)

class _CustomTargetHolder(ObjectHolder[_CT]):
    def __init__(self, target: _CT, interp: 'Interpreter'):
        super().__init__(target, interp)

    def __repr__(self) -> str:
        r = '<{} {}: {}>'
        h = self.held_object
        return r.format(self.__class__.__name__, h.get_id(), h.command)

    @noPosargs
    @noKwargs
    @InterpreterObject.method('full_path')
    def full_path_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self.interpreter.backend.get_target_filename_abs(self.held_object)

    @FeatureNew('custom_target.to_list', '0.54.0')
    @noPosargs
    @noKwargs
    @InterpreterObject.method('to_list')
    def to_list_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> T.List[build.CustomTargetIndex]:
        result = []
        for i in self.held_object:
            result.append(i)
        return result

    @noKwargs
    @typed_operator(MesonOperator.INDEX, int)
    @InterpreterObject.operator(MesonOperator.INDEX)
    def op_index(self, other: int) -> build.CustomTargetIndex:
        try:
            return self.held_object[other]
        except IndexError:
            raise InvalidArguments(f'Index {other} out of bounds of custom target {self.held_object.name} output of size {len(self.held_object)}.')

class CustomTargetHolder(_CustomTargetHolder[build.CustomTarget]):
    pass

class RunTargetHolder(ObjectHolder[build.RunTarget]):
    pass

class AliasTargetHolder(ObjectHolder[build.AliasTarget]):
    pass

class GeneratedListHolder(ObjectHolder[build.GeneratedList]):
    pass

class GeneratorHolder(ObjectHolder[build.Generator]):
    def __init__(self, gen: build.Generator, interpreter: 'Interpreter'):
        super().__init__(gen, interpreter)

    @typed_pos_args('generator.process', min_varargs=1, varargs=(str, mesonlib.File, build.CustomTarget, build.CustomTargetIndex, build.GeneratedList))
    @typed_kwargs(
        'generator.process',
        KwargInfo('preserve_path_from', (str, NoneType), since='0.45.0'),
        KwargInfo('extra_args', ContainerTypeInfo(list, str), listify=True, default=[]),
        ENV_KW.evolve(since='1.3.0')
    )
    @InterpreterObject.method('process')
    def process_method(self,
                       args: T.Tuple[T.List[T.Union[str, mesonlib.File, 'build.GeneratedTypes']]],
                       kwargs: 'kwargs.GeneratorProcess') -> build.GeneratedList:
        preserve_path_from = kwargs['preserve_path_from']
        if preserve_path_from is not None:
            preserve_path_from = os.path.normpath(preserve_path_from)
            if not os.path.isabs(preserve_path_from):
                # This is a bit of a hack. Fix properly before merging.
                raise InvalidArguments('Preserve_path_from must be an absolute path for now. Sorry.')

        if any(isinstance(a, (build.CustomTarget, build.CustomTargetIndex, build.GeneratedList)) for a in args[0]):
            FeatureNew.single_use(
                'Calling generator.process with CustomTarget or Index of CustomTarget.',
                '0.57.0', self.interpreter.subproject)

        gl = self.held_object.process_files(args[0], self.interpreter,
                                            preserve_path_from, extra_args=kwargs['extra_args'], env=kwargs['env'])

        return gl


class StructuredSourcesHolder(ObjectHolder[build.StructuredSources]):

    def __init__(self, sources: build.StructuredSources, interp: 'Interpreter'):
        super().__init__(sources, interp)

class OverrideExecutableHolder(BuildTargetHolder[build.OverrideExecutable]):
    @noPosargs
    @noKwargs
    @FeatureNew('OverrideExecutable.version', '1.9.0')
    @InterpreterObject.method('version')
    def version_method(self, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> str:
        return self.held_object.get_version(self.interpreter)
