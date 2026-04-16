# SPDX-License-Identifier: Apache-2.0
# Copyright 2013-2021 The Meson development team

from __future__ import annotations

import collections, importlib
import enum
import typing as T

from .base import DependencyCandidate, ExternalDependency, DependencyException, DependencyMethods, NotFoundDependency

from ..mesonlib import listify, PerMachine, MesonBugException, MesonException
from .. import mlog

if T.TYPE_CHECKING:
    from ..environment import Environment
    from .factory import DependencyFactory, DependencyGenerator, WrappedFactoryFunc
    from .base import DependencyObjectKWs

    TV_DepIDEntry = T.Union[str, bool, int, None, T.Tuple[str, ...]]
    TV_DepID = T.Tuple[T.Tuple[str, TV_DepIDEntry], ...]
    PackageTypes = T.Union[T.Type[ExternalDependency], DependencyFactory, DependencyCandidate, WrappedFactoryFunc]
    # Workaround for older python
    DependencyPackagesType = collections.UserDict[str, PackageTypes]
else:
    DependencyPackagesType = collections.UserDict

class DependencyPackages(DependencyPackagesType):
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

def get_dep_identifier(name: str, kwargs: DependencyObjectKWs) -> 'TV_DepID':
    identifier: 'TV_DepID' = (('name', name), )
    from ..interpreter.type_checking import DEPENDENCY_KWS
    nkwargs = T.cast('DependencyObjectKWs', {k.name: k.default for k in DEPENDENCY_KWS})
    nkwargs.update(kwargs)

    assert len(DEPENDENCY_KWS) == 20, \
           'Extra kwargs have been added to dependency(), please review if it makes sense to handle it here'
    for key, value in nkwargs.items():
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
                assert isinstance(i, str), i
            value = tuple(frozenset(listify(value)))
        elif isinstance(value, enum.Enum):
            value = value.value
            assert isinstance(value, str), 'for mypy'
        else:
            assert value is None or isinstance(value, (str, bool, int)), value
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

def find_external_dependency(name: str, env: 'Environment', kwargs: DependencyObjectKWs, candidates: T.Optional[T.List['DependencyGenerator']] = None) -> T.Union['ExternalDependency', NotFoundDependency]:
    assert name
    required = kwargs.get('required', True)
    lname = name.lower()
    if lname not in _packages_accept_language and kwargs.get('language') is not None:
        raise DependencyException(f'{name} dependency does not accept "language" keyword argument')

    # display the dependency name with correct casing
    display_name = display_name_map.get(lname, lname)

    for_machine = kwargs['native']
    type_text = PerMachine('Build-time', 'Run-time')[for_machine] + ' dependency'

    # build a list of dependency methods to try
    if candidates is None:
        candidates = _build_external_dependency_list(name, env, kwargs)

    pkg_exc: T.List[DependencyException] = []
    pkgdep:  T.List[ExternalDependency] = []
    details = ''
    tried_methods: T.List[str] = []

    for c in candidates:
        # try this dependency method
        try:
            d = c()
            d._check_version()
            pkgdep.append(d)
        except DependencyException as e:
            bettermsg = f'Dependency lookup for {name} with method {c.method!r} failed: {e}'
            mlog.debug(bettermsg)
            e.args = (bettermsg,)
            pkg_exc.append(e)
        except MesonException:
            raise
        except Exception as e:
            bettermsg = f'Dependency lookup for {name} with method {c.method!r} failed: {e}'
            raise MesonBugException(bettermsg) from e
        else:
            pkg_exc.append(None)
            details = d.log_details()
            if details:
                details = '(' + details + ') '
            if kwargs.get('language') is not None:
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
            tried_methods.append(c.method)

    # otherwise, the dependency could not be found
    tried = ' (tried {})'.format(mlog.format_list(tried_methods)) if tried_methods else ''
    mlog.log(type_text, mlog.bold(display_name), details + 'found:', mlog.red('NO'), tried)

    if required:
        # if an exception occurred with the first detection method, re-raise it
        # (on the grounds that it came from the preferred dependency detection
        # method)
        if pkg_exc and pkg_exc[0]:
            raise pkg_exc[0]

        # we have a list of failed ExternalDependency objects, so we can report
        # the methods we tried to find the dependency
        raise DependencyException(f'Dependency "{name}" not found' + tried)

    return NotFoundDependency(name, env)


def _build_external_dependency_list(name: str, env: 'Environment', kwargs: DependencyObjectKWs
                                    ) -> T.List['DependencyGenerator']:
    # Is there a specific dependency detector for this dependency?
    lname = name.lower()
    if lname in packages:
        entry = packages[lname]
        if isinstance(entry, type):
            if issubclass(entry, ExternalDependency):
                dep = [DependencyCandidate.from_dependency(name, entry, (env, kwargs))]
            else:
                raise MesonBugException(f'Got an invalid type in the dependency list: {entry!r}')
        elif isinstance(entry, DependencyCandidate):
            entry.arguments = (env, kwargs)
            dep = [entry]
        else:
            dep = entry(env, kwargs)
        return dep

    candidates: T.List['DependencyGenerator'] = []

    method = kwargs.get('method', DependencyMethods.AUTO)
    if method is DependencyMethods.AUTO:
        # Just use the standard detection methods.
        methods = [DependencyMethods.PKGCONFIG, DependencyMethods.EXTRAFRAMEWORK, DependencyMethods.CMAKE]
    else:
        # If it's explicitly requested, use that detection method (only).
        methods = [method]

    # Exclusive to when it is explicitly requested
    if DependencyMethods.DUB in methods:
        from .dub import DubDependency
        candidates.append(DependencyCandidate.from_dependency(name, DubDependency, (env, kwargs)))

    # Preferred first candidate for auto.
    if DependencyMethods.PKGCONFIG in methods:
        from .pkgconfig import PkgConfigDependency
        candidates.append(DependencyCandidate.from_dependency(name, PkgConfigDependency, (env, kwargs)))

    # On OSX only, try framework dependency detector.
    if DependencyMethods.EXTRAFRAMEWORK in methods:
        if env.machines[kwargs['native']].is_darwin():
            from .framework import ExtraFrameworkDependency
            candidates.append(DependencyCandidate.from_dependency(name, ExtraFrameworkDependency, (env, kwargs)))

    # Only use CMake:
    # - if it's explicitly requested
    # - as a last resort, since it might not work 100% (see #6113)
    if DependencyMethods.CMAKE in methods:
        from .cmake import CMakeDependency
        candidates.append(DependencyCandidate.from_dependency(name, CMakeDependency, (env, kwargs)))

    return candidates
