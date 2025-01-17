from __future__ import annotations

from .interpreterobjects import extract_required_kwarg
from .. import mlog
from .. import dependencies
from .. import build
from ..wrap import WrapMode
from ..mesonlib import extract_as_list, stringlistify, version_compare_many, listify
from ..options import OptionKey
from ..dependencies import Dependency, DependencyException, NotFoundDependency
from ..interpreterbase import (MesonInterpreterObject, FeatureNew,
                               InterpreterException, InvalidArguments)

import typing as T
if T.TYPE_CHECKING:
    from .interpreter import Interpreter
    from ..interpreterbase import TYPE_nkwargs, TYPE_nvar
    from .interpreterobjects import SubprojectHolder


class DependencyFallbacksHolder(MesonInterpreterObject):
    def __init__(self, interpreter: 'Interpreter', names: T.List[str], allow_fallback: T.Optional[bool] = None,
                 default_options: T.Optional[T.Dict[OptionKey, str]] = None) -> None:
        super().__init__(subproject=interpreter.subproject)
        self.interpreter = interpreter
        self.subproject = interpreter.subproject
        self.coredata = interpreter.coredata
        self.build = interpreter.build
        self.environment = interpreter.environment
        self.wrap_resolver = interpreter.environment.wrap_resolver
        self.allow_fallback = allow_fallback
        self.subproject_name: T.Optional[str] = None
        self.subproject_varname: T.Optional[str] = None
        self.subproject_kwargs = {'default_options': default_options or {}}
        self.names: T.List[str] = []
        self.forcefallback: bool = False
        self.nofallback: bool = False
        for name in names:
            if not name:
                raise InterpreterException('dependency_fallbacks empty name \'\' is not allowed')
            if '<' in name or '>' in name or '=' in name:
                raise InvalidArguments('Characters <, > and = are forbidden in dependency names. To specify'
                                       'version\n requirements use the \'version\' keyword argument instead.')
            if name in self.names:
                raise InterpreterException(f'dependency_fallbacks name {name!r} is duplicated')
            self.names.append(name)
        self._display_name = self.names[0] if self.names else '(anonymous)'

    def set_fallback(self, fbinfo: T.Optional[T.Union[T.List[str], str]]) -> None:
        # Legacy: This converts dependency()'s fallback kwargs.
        if fbinfo is None:
            return
        if self.allow_fallback is not None:
            raise InvalidArguments('"fallback" and "allow_fallback" arguments are mutually exclusive')
        fbinfo = stringlistify(fbinfo)
        if len(fbinfo) == 0:
            # dependency('foo', fallback: []) is the same as dependency('foo', allow_fallback: false)
            self.allow_fallback = False
            return
        if len(fbinfo) == 1:
            FeatureNew.single_use('Fallback without variable name', '0.53.0', self.subproject)
            subp_name, varname = fbinfo[0], None
        elif len(fbinfo) == 2:
            subp_name, varname = fbinfo
        else:
            raise InterpreterException('Fallback info must have one or two items.')
        self._subproject_impl(subp_name, varname)

    def _subproject_impl(self, subp_name: str, varname: str) -> None:
        assert self.subproject_name is None
        self.subproject_name = subp_name
        self.subproject_varname = varname

    def _do_dependency_cache(self, kwargs: TYPE_nkwargs, func_args: TYPE_nvar, func_kwargs: TYPE_nkwargs) -> T.Optional[Dependency]:
        name = func_args[0]
        cached_dep = self._get_cached_dep(name, kwargs)
        if cached_dep:
            self._verify_fallback_consistency(cached_dep)
        return cached_dep

    def _do_dependency(self, kwargs: TYPE_nkwargs, func_args: TYPE_nvar, func_kwargs: TYPE_nkwargs) -> T.Optional[Dependency]:
        # Note that there is no df.dependency() method, this is called for names
        # given as positional arguments to dependency_fallbacks(name1, ...).
        # We use kwargs from the dependency() function, for things like version,
        # module, etc.
        name = func_args[0]
        self._handle_featurenew_dependencies(name)
        dep = dependencies.find_external_dependency(name, self.environment, kwargs)
        if dep.found():
            for_machine = self.interpreter.machine_from_native_kwarg(kwargs)
            identifier = dependencies.get_dep_identifier(name, kwargs)
            self.coredata.deps[for_machine].put(identifier, dep)
            return dep
        return None

    def _do_existing_subproject(self, kwargs: TYPE_nkwargs, func_args: TYPE_nvar, func_kwargs: TYPE_nkwargs) -> T.Optional[Dependency]:
        subp_name = func_args[0]
        varname = self.subproject_varname
        if subp_name and self._get_subproject(subp_name):
            return self._get_subproject_dep(subp_name, varname, kwargs)
        return None

    def _do_subproject(self, kwargs: TYPE_nkwargs, func_args: TYPE_nvar, func_kwargs: TYPE_nkwargs) -> T.Optional[Dependency]:
        if self.forcefallback:
            mlog.log('Looking for a fallback subproject for the dependency',
                     mlog.bold(self._display_name), 'because:\nUse of fallback dependencies is forced.')
        elif self.nofallback:
            mlog.log('Not looking for a fallback subproject for the dependency',
                     mlog.bold(self._display_name), 'because:\nUse of fallback dependencies is disabled.')
            return None
        else:
            mlog.log('Looking for a fallback subproject for the dependency',
                     mlog.bold(self._display_name))

        # dependency('foo', static: true) should implicitly add
        # default_options: ['default_library=static']
        static = kwargs.get('static')
        default_options = func_kwargs.get('default_options', {})
        if static is not None and 'default_library' not in default_options:
            default_library = 'static' if static else 'shared'
            mlog.log(f'Building fallback subproject with default_library={default_library}')
            default_options[OptionKey('default_library')] = default_library
            func_kwargs['default_options'] = default_options

        # Configure the subproject
        subp_name = self.subproject_name
        varname = self.subproject_varname
        func_kwargs.setdefault('version', [])
        if 'default_options' in kwargs and isinstance(kwargs['default_options'], str):
            func_kwargs['default_options'] = listify(kwargs['default_options'])
        self.interpreter.do_subproject(subp_name, func_kwargs)
        return self._get_subproject_dep(subp_name, varname, kwargs)

    def _get_subproject(self, subp_name: str) -> T.Optional[SubprojectHolder]:
        sub = self.interpreter.subprojects.get(subp_name)
        if sub and sub.found():
            return sub
        return None

    def _get_subproject_dep(self, subp_name: str, varname: str, kwargs: TYPE_nkwargs) -> T.Optional[Dependency]:
        # Verify the subproject is found
        subproject = self._get_subproject(subp_name)
        if not subproject:
            mlog.log('Dependency', mlog.bold(self._display_name), 'from subproject',
                     mlog.bold(subp_name), 'found:', mlog.red('NO'),
                     mlog.blue('(subproject failed to configure)'))
            return None

        # The subproject has been configured. If for any reason the dependency
        # cannot be found in this subproject we have to return not-found object
        # instead of None, because we don't want to continue the lookup on the
        # system.

        # Check if the subproject overridden at least one of the names we got.
        cached_dep = None
        for name in self.names:
            cached_dep = self._get_cached_dep(name, kwargs)
            if cached_dep:
                break

        # If we have cached_dep we did all the checks and logging already in
        # self._get_cached_dep().
        if cached_dep:
            self._verify_fallback_consistency(cached_dep)
            return cached_dep

        # Legacy: Use the variable name if provided instead of relying on the
        # subproject to override one of our dependency names
        if not varname:
            # If no variable name is specified, check if the wrap file has one.
            # If the wrap file has a variable name, better use it because the
            # subproject most probably is not using meson.override_dependency().
            for name in self.names:
                varname = self.wrap_resolver.get_varname(subp_name, name)
                if varname:
                    break
        if not varname:
            mlog.warning(f'Subproject {subp_name!r} did not override {self._display_name!r} dependency and no variable name specified')
            mlog.log('Dependency', mlog.bold(self._display_name), 'from subproject',
                     mlog.bold(subproject.subdir), 'found:', mlog.red('NO'))
            return self._notfound_dependency()

        var_dep = self._get_subproject_variable(subproject, varname) or self._notfound_dependency()
        if not var_dep.found():
            mlog.log('Dependency', mlog.bold(self._display_name), 'from subproject',
                     mlog.bold(subproject.subdir), 'found:', mlog.red('NO'))
            return var_dep

        wanted = stringlistify(kwargs.get('version', []))
        found = var_dep.get_version()
        if not self._check_version(wanted, found):
            mlog.log('Dependency', mlog.bold(self._display_name), 'from subproject',
                     mlog.bold(subproject.subdir), 'found:', mlog.red('NO'),
                     'found', mlog.normal_cyan(found), 'but need:',
                     mlog.bold(', '.join([f"'{e}'" for e in wanted])))
            return self._notfound_dependency()

        mlog.log('Dependency', mlog.bold(self._display_name), 'from subproject',
                 mlog.bold(subproject.subdir), 'found:', mlog.green('YES'),
                 mlog.normal_cyan(found) if found else None)
        return var_dep

    def _get_cached_dep(self, name: str, kwargs: TYPE_nkwargs) -> T.Optional[Dependency]:
        # Unlike other methods, this one returns not-found dependency instead
        # of None in the case the dependency is cached as not-found, or if cached
        # version does not match. In that case we don't want to continue with
        # other candidates.
        for_machine = self.interpreter.machine_from_native_kwarg(kwargs)
        identifier = dependencies.get_dep_identifier(name, kwargs)
        wanted_vers = stringlistify(kwargs.get('version', []))

        override = self.build.dependency_overrides[for_machine].get(identifier)
        if override:
            info = [mlog.blue('(overridden)' if override.explicit else '(cached)')]
            cached_dep = override.dep
            # We don't implicitly override not-found dependencies, but user could
            # have explicitly called meson.override_dependency() with a not-found
            # dep.
            if not cached_dep.found():
                mlog.log('Dependency', mlog.bold(self._display_name),
                         'found:', mlog.red('NO'), *info)
                return cached_dep
        elif self.forcefallback and self.subproject_name:
            cached_dep = None
        else:
            info = [mlog.blue('(cached)')]
            cached_dep = self.coredata.deps[for_machine].get(identifier)

        if cached_dep:
            found_vers = cached_dep.get_version()
            if not self._check_version(wanted_vers, found_vers):
                if not override:
                    # We cached this dependency on disk from a previous run,
                    # but it could got updated on the system in the meantime.
                    return None
                mlog.log('Dependency', mlog.bold(name),
                         'found:', mlog.red('NO'),
                         'found', mlog.normal_cyan(found_vers), 'but need:',
                         mlog.bold(', '.join([f"'{e}'" for e in wanted_vers])),
                         *info)
                return self._notfound_dependency()
            if found_vers:
                info = [mlog.normal_cyan(found_vers), *info]
            mlog.log('Dependency', mlog.bold(self._display_name),
                     'found:', mlog.green('YES'), *info)
            return cached_dep
        return None

    def _get_subproject_variable(self, subproject: SubprojectHolder, varname: str) -> T.Optional[Dependency]:
        try:
            var_dep = subproject.get_variable_method([varname], {})
        except InvalidArguments:
            var_dep = None
        if not isinstance(var_dep, Dependency):
            mlog.warning(f'Variable {varname!r} in the subproject {subproject.subdir!r} is',
                         'not found' if var_dep is None else 'not a dependency object')
            return None
        return var_dep

    def _verify_fallback_consistency(self, cached_dep: Dependency) -> None:
        subp_name = self.subproject_name
        varname = self.subproject_varname
        subproject = self._get_subproject(subp_name)
        if subproject and varname:
            var_dep = self._get_subproject_variable(subproject, varname)
            if var_dep and cached_dep.found() and var_dep != cached_dep:
                mlog.warning(f'Inconsistency: Subproject has overridden the dependency with another variable than {varname!r}')

    def _handle_featurenew_dependencies(self, name: str) -> None:
        'Do a feature check on dependencies used by this subproject'
        if name == 'mpi':
            FeatureNew.single_use('MPI Dependency', '0.42.0', self.subproject)
        elif name == 'pcap':
            FeatureNew.single_use('Pcap Dependency', '0.42.0', self.subproject)
        elif name == 'vulkan':
            FeatureNew.single_use('Vulkan Dependency', '0.42.0', self.subproject)
        elif name == 'libwmf':
            FeatureNew.single_use('LibWMF Dependency', '0.44.0', self.subproject)
        elif name == 'openmp':
            FeatureNew.single_use('OpenMP Dependency', '0.46.0', self.subproject)

    def _notfound_dependency(self) -> NotFoundDependency:
        return NotFoundDependency(self.names[0] if self.names else '', self.environment)

    @staticmethod
    def _check_version(wanted: T.List[str], found: str) -> bool:
        if not wanted:
            return True
        return not (found == 'undefined' or not version_compare_many(found, wanted)[0])

    def _get_candidates(self) -> T.List[T.Tuple[T.Callable[[TYPE_nkwargs, TYPE_nvar, TYPE_nkwargs], T.Optional[Dependency]], TYPE_nvar, TYPE_nkwargs]]:
        candidates = []
        # 1. check if any of the names is cached already.
        for name in self.names:
            candidates.append((self._do_dependency_cache, [name], {}))
        # 2. check if the subproject fallback has already been configured.
        if self.subproject_name:
            candidates.append((self._do_existing_subproject, [self.subproject_name], self.subproject_kwargs))
        # 3. check external dependency if we are not forced to use subproject
        if not self.forcefallback or not self.subproject_name:
            for name in self.names:
                candidates.append((self._do_dependency, [name], {}))
        # 4. configure the subproject
        if self.subproject_name:
            candidates.append((self._do_subproject, [self.subproject_name], self.subproject_kwargs))
        return candidates

    def lookup(self, kwargs: TYPE_nkwargs, force_fallback: bool = False) -> Dependency:
        mods = extract_as_list(kwargs, 'modules')
        if mods:
            self._display_name += ' (modules: {})'.format(', '.join(str(i) for i in mods))

        disabled, required, feature = extract_required_kwarg(kwargs, self.subproject)
        if disabled:
            mlog.log('Dependency', mlog.bold(self._display_name), 'skipped: feature', mlog.bold(feature), 'disabled')
            return self._notfound_dependency()

        # Check if usage of the subproject fallback is forced
        wrap_mode = WrapMode.from_string(self.coredata.get_option(OptionKey('wrap_mode')))
        force_fallback_for = self.coredata.get_option(OptionKey('force_fallback_for'))
        assert isinstance(force_fallback_for, list), 'for mypy'
        self.nofallback = wrap_mode == WrapMode.nofallback
        self.forcefallback = (force_fallback or
                              wrap_mode == WrapMode.forcefallback or
                              any(name in force_fallback_for for name in self.names) or
                              self.subproject_name in force_fallback_for)

        # Add an implicit subproject fallback if none has been set explicitly,
        # unless implicit fallback is not allowed.
        # Legacy: self.allow_fallback can be None when that kwarg is not defined
        # in dependency('name'). In that case we don't want to use implicit
        # fallback when required is false because user will typically fallback
        # manually using cc.find_library() for example.
        if not self.subproject_name and self.allow_fallback is not False:
            for name in self.names:
                subp_name, varname = self.wrap_resolver.find_dep_provider(name)
                if subp_name:
                    self.forcefallback |= subp_name in force_fallback_for
                    if self.forcefallback or self.allow_fallback is True or required or self._get_subproject(subp_name):
                        self._subproject_impl(subp_name, varname)
                    break

        candidates = self._get_candidates()

        # writing just "dependency('')" is an error, because it can only fail
        if not candidates and required:
            raise InvalidArguments('Dependency is required but has no candidates.')

        # Try all candidates, only the last one is really required.
        last = len(candidates) - 1
        for i, item in enumerate(candidates):
            func, func_args, func_kwargs = item
            func_kwargs['required'] = required and (i == last)
            kwargs['required'] = required and (i == last)
            dep = func(kwargs, func_args, func_kwargs)
            if dep and dep.found():
                # Override this dependency to have consistent results in subsequent
                # dependency lookups.
                for name in self.names:
                    for_machine = self.interpreter.machine_from_native_kwarg(kwargs)
                    identifier = dependencies.get_dep_identifier(name, kwargs)
                    if identifier not in self.build.dependency_overrides[for_machine]:
                        self.build.dependency_overrides[for_machine][identifier] = \
                            build.DependencyOverride(dep, self.interpreter.current_node, explicit=False)
                return dep
            elif required and (dep or i == last):
                # This was the last candidate or the dependency has been cached
                # as not-found, or cached dependency version does not match,
                # otherwise func() would have returned None instead.
                raise DependencyException(f'Dependency {self._display_name!r} is required but not found.')
            elif dep:
                # Same as above, but the dependency is not required.
                return dep
        return self._notfound_dependency()
