# SPDX-License-Identifier: Apache-2.0
# Copyright 2013-2021 The Meson development team

from __future__ import annotations

import collections, functools, importlib
import typing as T

from .base import ExternalDependency, DependencyException, DependencyMethods, NotFoundDependency

from ..mesonlib import listify, MachineChoice, PerMachine
from .. import mlog

if T.TYPE_CHECKING:
    from ..environment import Environment
    from .factory import DependencyFactory, WrappedFactoryFunc, DependencyGenerator

    TV_DepIDEntry = T.Union[str, bool, int, T.Tuple[str, ...]]
    TV_DepID = T.Tuple[T.Tuple[str, TV_DepIDEntry], ...]
    PackageTypes = T.Union[T.Type[ExternalDependency], DependencyFactory, WrappedFactoryFunc]

class DependencyPackages(collections.UserDict):
    data: T.Dict[str, PackageTypes]
    defaults: T.Dict[str, str] = {}

    def __missing__(self, key: str) -> PackageTypes:
        if key in self.defaults:
            modn = self.defaults[key]
            importlib.import_module(f'mesonbuild.dependencies.{modn}')

            return self.data[key]
        raise KeyError(key)

    def __contains__(self, key: object) -> bool:
        return key in self.defaults or key in self.data

# These must be defined in this file to avoid cyclical references.
packages = DependencyPackages()
_packages_accept_language: T.Set[str] = set()

def get_dep_identifier(name: str, kwargs: T.Dict[str, T.Any]) -> 'TV_DepID':
    identifier: 'TV_DepID' = (('name', name), )
    from ..interpreter import permitted_dependency_kwargs
    assert len(permitted_dependency_kwargs) == 19, \
           'Extra kwargs have been added to dependency(), please review if it makes sense to handle it here'
    for key, value in kwargs.items():
        # 'version' is irrelevant for caching; the caller must check version matches
        # 'native' is handled above with `for_machine`
        # 'required' is irrelevant for caching; the caller handles it separately
        # 'fallback' and 'allow_fallback' is not part of the cache because,
        #     once a dependency has been found through a fallback, it should
        #     be used for the rest of the Meson run.
        # 'default_options' is only used in fallback case
        # 'not_found_message' has no impact on the dependency lookup
        # 'include_type' is handled after the dependency lookup
        if key in {'version', 'native', 'required', 'fallback', 'allow_fallback', 'default_options',
                   'not_found_message', 'include_type'}:
            continue
        # All keyword arguments are strings, ints, or lists (or lists of lists)
        if isinstance(value, list):
            for i in value:
                assert isinstance(i, str)
            value = tuple(frozenset(listify(value)))
        else:
            assert isinstance(value, (str, bool, int))
        identifier = (*identifier, (key, value),)
    return identifier

display_name_map = {
    'boost': 'Boost',
    'cuda': 'CUDA',
    'dub': 'DUB',
    'gmock': 'GMock',
    'gtest': 'GTest',
    'hdf5': 'HDF5',
    'llvm': 'LLVM',
    'mpi': 'MPI',
    'netcdf': 'NetCDF',
    'openmp': 'OpenMP',
    'wxwidgets': 'WxWidgets',
}

def find_external_dependency(name: str, env: 'Environment', kwargs: T.Dict[str, object], candidates: T.Optional[T.List['DependencyGenerator']] = None) -> T.Union['ExternalDependency', NotFoundDependency]:
    assert name
    required = kwargs.get('required', True)
    if not isinstance(required, bool):
        raise DependencyException('Keyword "required" must be a boolean.')
    if not isinstance(kwargs.get('method', ''), str):
        raise DependencyException('Keyword "method" must be a string.')
    lname = name.lower()
    if lname not in _packages_accept_language and 'language' in kwargs:
        raise DependencyException(f'{name} dependency does not accept "language" keyword argument')
    if not isinstance(kwargs.get('version', ''), (str, list)):
        raise DependencyException('Keyword "Version" must be string or list.')

    # display the dependency name with correct casing
    display_name = display_name_map.get(lname, lname)

    for_machine = MachineChoice.BUILD if kwargs.get('native', False) else MachineChoice.HOST

    type_text = PerMachine('Build-time', 'Run-time')[for_machine] + ' dependency'

    # build a list of dependency methods to try
    if candidates is None:
        candidates = _build_external_dependency_list(name, env, for_machine, kwargs)

    pkg_exc: T.List[DependencyException] = []
    pkgdep:  T.List[ExternalDependency] = []
    details = ''

    for c in candidates:
        # try this dependency method
        try:
            d = c()
            d._check_version()
            pkgdep.append(d)
        except DependencyException as e:
            assert isinstance(c, functools.partial), 'for mypy'
            bettermsg = f'Dependency lookup for {name} with method {c.func.log_tried()!r} failed: {e}'
            mlog.debug(bettermsg)
            e.args = (bettermsg,)
            pkg_exc.append(e)
        else:
            pkg_exc.append(None)
            details = d.log_details()
            if details:
                details = '(' + details + ') '
            if 'language' in kwargs:
                details += 'for ' + d.language + ' '

            # if the dependency was found
            if d.found():

                info: mlog.TV_LoggableList = []
                if d.version:
                    info.append(mlog.normal_cyan(d.version))

                log_info = d.log_info()
                if log_info:
                    info.append('(' + log_info + ')')

                mlog.log(type_text, mlog.bold(display_name), details + 'found:', mlog.green('YES'), *info)

                return d

    # otherwise, the dependency could not be found
    tried_methods = [d.log_tried() for d in pkgdep if d.log_tried()]
    if tried_methods:
        tried = mlog.format_list(tried_methods)
    else:
        tried = ''

    mlog.log(type_text, mlog.bold(display_name), details + 'found:', mlog.red('NO'),
             f'(tried {tried})' if tried else '')

    if required:
        # if an exception occurred with the first detection method, re-raise it
        # (on the grounds that it came from the preferred dependency detection
        # method)
        if pkg_exc and pkg_exc[0]:
            raise pkg_exc[0]

        # we have a list of failed ExternalDependency objects, so we can report
        # the methods we tried to find the dependency
        raise DependencyException(f'Dependency "{name}" not found' +
                                  (f', tried {tried}' if tried else ''))

    return NotFoundDependency(name, env)


def _build_external_dependency_list(name: str, env: 'Environment', for_machine: MachineChoice,
                                    kwargs: T.Dict[str, T.Any]) -> T.List['DependencyGenerator']:
    # First check if the method is valid
    if 'method' in kwargs and kwargs['method'] not in [e.value for e in DependencyMethods]:
        raise DependencyException('method {!r} is invalid'.format(kwargs['method']))

    # Is there a specific dependency detector for this dependency?
    lname = name.lower()
    if lname in packages:
        # Create the list of dependency object constructors using a factory
        # class method, if one exists, otherwise the list just consists of the
        # constructor
        if isinstance(packages[lname], type):
            entry1 = T.cast('T.Type[ExternalDependency]', packages[lname])  # mypy doesn't understand isinstance(..., type)
            if issubclass(entry1, ExternalDependency):
                func: T.Callable[[], 'ExternalDependency'] = functools.partial(entry1, env, kwargs)
                dep = [func]
        else:
            entry2 = T.cast('T.Union[DependencyFactory, WrappedFactoryFunc]', packages[lname])
            dep = entry2(env, for_machine, kwargs)
        return dep

    candidates: T.List['DependencyGenerator'] = []

    if kwargs.get('method', 'auto') == 'auto':
        # Just use the standard detection methods.
        methods = ['pkg-config', 'extraframework', 'cmake']
    else:
        # If it's explicitly requested, use that detection method (only).
        methods = [kwargs['method']]

    # Exclusive to when it is explicitly requested
    if 'dub' in methods:
        from .dub import DubDependency
        candidates.append(functools.partial(DubDependency, name, env, kwargs))

    # Preferred first candidate for auto.
    if 'pkg-config' in methods:
        from .pkgconfig import PkgConfigDependency
        candidates.append(functools.partial(PkgConfigDependency, name, env, kwargs))

    # On OSX only, try framework dependency detector.
    if 'extraframework' in methods:
        if env.machines[for_machine].is_darwin():
            from .framework import ExtraFrameworkDependency
            candidates.append(functools.partial(ExtraFrameworkDependency, name, env, kwargs))

    # Only use CMake:
    # - if it's explicitly requested
    # - as a last resort, since it might not work 100% (see #6113)
    if 'cmake' in methods:
        from .cmake import CMakeDependency
        candidates.append(functools.partial(CMakeDependency, name, env, kwargs))

    return candidates
