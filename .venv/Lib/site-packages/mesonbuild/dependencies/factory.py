# SPDX-License-Identifier: Apache-2.0
# Copyright 2013-2021 The Meson development team
# Copyright © 2021-2025 Intel Corporation

from __future__ import annotations

import functools
import typing as T

from ..mesonlib import MachineChoice
from .base import DependencyCandidate, DependencyException, DependencyMethods
from .base import process_method_kw
from .base import BuiltinDependency, SystemDependency
from .cmake import CMakeDependency
from .framework import ExtraFrameworkDependency
from .pkgconfig import PkgConfigDependency

if T.TYPE_CHECKING:
    from typing_extensions import TypeAlias

    from .base import DependencyObjectKWs, ExternalDependency, DepType
    from .configtool import ConfigToolDependency
    from ..environment import Environment

    # TODO: remove this?
    DependencyGenerator: TypeAlias = DependencyCandidate[ExternalDependency]
    FactoryFunc = T.Callable[
        [
            'Environment',
            DependencyObjectKWs,
            T.List[DependencyMethods]
        ],
        T.List[DependencyGenerator]
    ]

    WrappedFactoryFunc = T.Callable[
        [
            'Environment',
            DependencyObjectKWs,
        ],
        T.List[DependencyGenerator]
    ]

class DependencyFactory:

    """Factory to get dependencies from multiple sources.

    This class provides an initializer that takes a set of names and classes
    for various kinds of dependencies. When the initialized object is called
    it returns a list of callables return Dependency objects to try in order.

    :param name: The name of the dependency. This will be passed as the name
        parameter of the each dependency unless it is overridden on a per
        type basis.
    :param methods: An ordered list of DependencyMethods. This is the order
        dependencies will be returned in unless they are removed by the
        _process_method function
    :param extra_kwargs: Additional keyword arguments to add when creating the
        DependencyCandidate
    :param pkgconfig: A custom PackageConfig lookup to use
    :param cmake: A custom CMake lookup to use
    :param framework: A custom AppleFramework lookup to use
    :param configtool: A custom ConfigTool lookup to use. If
        DependencyMethods.CONFIG_TOOL is in the `:param:methods` argument,
        this must be set.
    :param builtin: A custom Builtin lookup to use. If
        DependencyMethods.BUILTIN is in the `:param:methods` argument,
        this must be set.
    :param system: A custom System lookup to use. If
        DependencyMethods.SYSTEM is in the `:param:methods` argument,
        this must be set.
    """

    def __init__(self, name: str, methods: T.List[DependencyMethods], *,
                 extra_kwargs: T.Optional[DependencyObjectKWs] = None,
                 pkgconfig: T.Union[DependencyCandidate[PkgConfigDependency], T.Type[PkgConfigDependency], None] = PkgConfigDependency,
                 cmake: T.Union[DependencyCandidate[CMakeDependency], T.Type[CMakeDependency], None] = CMakeDependency,
                 framework: T.Union[DependencyCandidate[ExtraFrameworkDependency], T.Type[ExtraFrameworkDependency], None] = ExtraFrameworkDependency,
                 configtool: T.Union[DependencyCandidate[ConfigToolDependency], T.Type[ConfigToolDependency], None] = None,
                 builtin: T.Union[DependencyCandidate[BuiltinDependency], T.Type[BuiltinDependency], None] = None,
                 system: T.Union[DependencyCandidate[SystemDependency], T.Type[SystemDependency], None] = None):

        if DependencyMethods.CONFIG_TOOL in methods and not configtool:
            raise DependencyException('A configtool dependency must have a custom class')
        if DependencyMethods.BUILTIN in methods and not builtin:
            raise DependencyException('A builtin dependency must have a custom class')
        if DependencyMethods.SYSTEM in methods and not system:
            raise DependencyException('A system dependency must have a custom class')

        def make(arg: T.Union[DependencyCandidate[DepType], T.Type[DepType], None]) -> T.Optional[DependencyCandidate[DepType]]:
            if arg is None or isinstance(arg, DependencyCandidate):
                return arg
            return DependencyCandidate.from_dependency(name, arg)

        self.extra_kwargs = extra_kwargs
        self.methods = methods
        self.classes: T.Mapping[DependencyMethods, T.Optional[DependencyCandidate[ExternalDependency]]] = {
            # Just attach the correct name right now, either the generic name
            # or the method specific name.
            DependencyMethods.EXTRAFRAMEWORK: make(framework),
            DependencyMethods.PKGCONFIG: make(pkgconfig),
            DependencyMethods.CMAKE: make(cmake),
            DependencyMethods.SYSTEM: make(system),
            DependencyMethods.BUILTIN: make(builtin),
            DependencyMethods.CONFIG_TOOL: make(configtool),
        }

    @staticmethod
    def _process_method(method: DependencyMethods, env: 'Environment', for_machine: MachineChoice) -> bool:
        """Report whether a method is valid or not.

        If the method is valid, return true, otherwise return false. This is
        used in a list comprehension to filter methods that are not possible.

        By default this only remove EXTRAFRAMEWORK dependencies for non-mac platforms.
        """
        # Extra frameworks are only valid for macOS and other apple products
        if (method is DependencyMethods.EXTRAFRAMEWORK and
                not env.machines[for_machine].is_darwin()):
            return False
        return True

    def __call__(self, env: 'Environment', kwargs: DependencyObjectKWs) -> T.List['DependencyGenerator']:
        """Return a list of Dependencies with the arguments already attached."""
        methods = process_method_kw(self.methods, kwargs)
        if self.extra_kwargs:
            nwargs = self.extra_kwargs.copy()
            nwargs.update(kwargs)
        else:
            nwargs = kwargs.copy()

        ret: T.List[DependencyGenerator] = []
        for m in methods:
            if self._process_method(m, env, kwargs['native']):
                c = self.classes[m]
                if c is None:
                    continue
                c.arguments = (env, nwargs)
                ret.append(c)
        return ret


def factory_methods(methods: T.Set[DependencyMethods]) -> T.Callable[['FactoryFunc'], 'WrappedFactoryFunc']:
    """Decorator for handling methods for dependency factory functions.

    This helps to make factory functions self documenting
    >>> @factory_methods([DependencyMethods.PKGCONFIG, DependencyMethods.CMAKE])
    >>> def factory(env: Environment, for_machine: MachineChoice, kwargs: DependencyObjectKWs, methods: T.List[DependencyMethods]) -> T.List['DependencyGenerator']:
    >>>     pass
    """

    def inner(func: 'FactoryFunc') -> 'WrappedFactoryFunc':

        @functools.wraps(func)
        def wrapped(env: 'Environment', kwargs: DependencyObjectKWs) -> T.List['DependencyGenerator']:
            return func(env, kwargs, process_method_kw(methods, kwargs))

        return wrapped

    return inner
