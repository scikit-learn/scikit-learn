# SPDX-License-Identifier: Apache-2.0
# Copyright 2013-2021 The Meson development team
# Copyright © 2021-2023 Intel Corporation

from __future__ import annotations

import functools
import typing as T

from .base import DependencyException, DependencyMethods
from .base import process_method_kw
from .base import BuiltinDependency, SystemDependency
from .cmake import CMakeDependency
from .framework import ExtraFrameworkDependency
from .pkgconfig import PkgConfigDependency

if T.TYPE_CHECKING:
    from .base import ExternalDependency
    from .configtool import ConfigToolDependency
    from ..environment import Environment
    from ..mesonlib import MachineChoice

    DependencyGenerator = T.Callable[[], ExternalDependency]
    FactoryFunc = T.Callable[
        [
            'Environment',
            MachineChoice,
            T.Dict[str, T.Any],
            T.List[DependencyMethods]
        ],
        T.List[DependencyGenerator]
    ]

    WrappedFactoryFunc = T.Callable[
        [
            'Environment',
            MachineChoice,
            T.Dict[str, T.Any]
        ],
        T.List[DependencyGenerator]
    ]

    # This should be str, Environment, T.Dict[str, T.Any], T.Optional[str]
    # But if you try that, you get error: Cannot infer type of lambda
    CmakeDependencyFunc = T.Callable[..., CMakeDependency]

class DependencyFactory:

    """Factory to get dependencies from multiple sources.

    This class provides an initializer that takes a set of names and classes
    for various kinds of dependencies. When the initialized object is called
    it returns a list of callables return Dependency objects to try in order.

    :name: The name of the dependency. This will be passed as the name
        parameter of the each dependency unless it is overridden on a per
        type basis.
    :methods: An ordered list of DependencyMethods. This is the order
        dependencies will be returned in unless they are removed by the
        _process_method function
    :*_name: This will overwrite the name passed to the corresponding class.
        For example, if the name is 'zlib', but cmake calls the dependency
        'Z', then using `cmake_name='Z'` will pass the name as 'Z' to cmake.
    :*_class: A *type* or callable that creates a class, and has the
        signature of an ExternalDependency
    :system_class: If you pass DependencyMethods.SYSTEM in methods, you must
        set this argument.
    """

    def __init__(self, name: str, methods: T.List[DependencyMethods], *,
                 extra_kwargs: T.Optional[T.Dict[str, T.Any]] = None,
                 pkgconfig_name: T.Optional[str] = None,
                 pkgconfig_class: 'T.Type[PkgConfigDependency]' = PkgConfigDependency,
                 cmake_name: T.Optional[str] = None,
                 cmake_class: 'T.Union[T.Type[CMakeDependency], CmakeDependencyFunc]' = CMakeDependency,
                 configtool_class: 'T.Optional[T.Type[ConfigToolDependency]]' = None,
                 framework_name: T.Optional[str] = None,
                 framework_class: 'T.Type[ExtraFrameworkDependency]' = ExtraFrameworkDependency,
                 builtin_class: 'T.Type[BuiltinDependency]' = BuiltinDependency,
                 system_class: 'T.Type[SystemDependency]' = SystemDependency):

        if DependencyMethods.CONFIG_TOOL in methods and not configtool_class:
            raise DependencyException('A configtool must have a custom class')

        self.extra_kwargs = extra_kwargs or {}
        self.methods = methods
        self.classes: T.Dict[
            DependencyMethods,
            T.Callable[['Environment', T.Dict[str, T.Any]], ExternalDependency]
        ] = {
            # Just attach the correct name right now, either the generic name
            # or the method specific name.
            DependencyMethods.EXTRAFRAMEWORK: functools.partial(framework_class, framework_name or name),
            DependencyMethods.PKGCONFIG: functools.partial(pkgconfig_class, pkgconfig_name or name),
            DependencyMethods.CMAKE: functools.partial(cmake_class, cmake_name or name),
            DependencyMethods.SYSTEM: functools.partial(system_class, name),
            DependencyMethods.BUILTIN: functools.partial(builtin_class, name),
            DependencyMethods.CONFIG_TOOL: None,
        }
        if configtool_class is not None:
            self.classes[DependencyMethods.CONFIG_TOOL] = functools.partial(configtool_class, name)

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

    def __call__(self, env: 'Environment', for_machine: MachineChoice,
                 kwargs: T.Dict[str, T.Any]) -> T.List['DependencyGenerator']:
        """Return a list of Dependencies with the arguments already attached."""
        methods = process_method_kw(self.methods, kwargs)
        nwargs = self.extra_kwargs.copy()
        nwargs.update(kwargs)

        return [functools.partial(self.classes[m], env, nwargs) for m in methods
                if self._process_method(m, env, for_machine)]


def factory_methods(methods: T.Set[DependencyMethods]) -> T.Callable[['FactoryFunc'], 'WrappedFactoryFunc']:
    """Decorator for handling methods for dependency factory functions.

    This helps to make factory functions self documenting
    >>> @factory_methods([DependencyMethods.PKGCONFIG, DependencyMethods.CMAKE])
    >>> def factory(env: Environment, for_machine: MachineChoice, kwargs: T.Dict[str, T.Any], methods: T.List[DependencyMethods]) -> T.List['DependencyGenerator']:
    >>>     pass
    """

    def inner(func: 'FactoryFunc') -> 'WrappedFactoryFunc':

        @functools.wraps(func)
        def wrapped(env: 'Environment', for_machine: MachineChoice, kwargs: T.Dict[str, T.Any]) -> T.List['DependencyGenerator']:
            return func(env, for_machine, kwargs, process_method_kw(methods, kwargs))

        return wrapped

    return inner
