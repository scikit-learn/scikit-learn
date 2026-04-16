# SPDX-License-Identifier: Apache-2.0
# Copyright 2013-2019 The Meson development team

from __future__ import annotations

import typing as T

from .base import DependencyCandidate, DependencyMethods, detect_compiler, SystemDependency
from .cmake import CMakeDependency
from .detect import packages
from .pkgconfig import PkgConfigDependency
from .factory import factory_methods

if T.TYPE_CHECKING:
    from . factory import DependencyGenerator
    from ..environment import Environment
    from .base import DependencyObjectKWs


@factory_methods({DependencyMethods.PKGCONFIG, DependencyMethods.CMAKE, DependencyMethods.SYSTEM})
def coarray_factory(env: 'Environment',
                    kwargs: DependencyObjectKWs,
                    methods: T.List[DependencyMethods]) -> T.List['DependencyGenerator']:
    kwargs['language'] = 'fortran'
    for_machine = kwargs['native']
    fcid = detect_compiler('coarray', env, for_machine, 'fortran').get_id()
    candidates: T.List['DependencyGenerator'] = []

    if fcid == 'gcc':
        # OpenCoarrays is the most commonly used method for Fortran Coarray with GCC
        if DependencyMethods.PKGCONFIG in methods:
            for pkg in ['caf-openmpi', 'caf']:
                candidates.append(DependencyCandidate.from_dependency(
                    pkg, PkgConfigDependency, (env, kwargs)))

        if DependencyMethods.CMAKE in methods:
            nkwargs = kwargs
            if not kwargs.get('modules'):
                nkwargs = kwargs.copy()
                nkwargs['modules'] = ['OpenCoarrays::caf_mpi']
            candidates.append(DependencyCandidate.from_dependency(
                'OpenCoarrays', CMakeDependency, (env, nkwargs)))

    if DependencyMethods.SYSTEM in methods:
        candidates.append(DependencyCandidate.from_dependency(
            'coarray', CoarrayDependency, (env, kwargs)))

    return candidates


packages['coarray'] = coarray_factory


class CoarrayDependency(SystemDependency):
    """
    Coarrays are a Fortran 2008 feature.

    Coarrays are sometimes implemented via external library (GCC+OpenCoarrays),
    while other compilers just build in support (Cray, IBM, Intel, NAG).
    Coarrays may be thought of as a high-level language abstraction of
    low-level MPI calls.
    """
    def __init__(self, name: str, environment: 'Environment', kwargs: DependencyObjectKWs) -> None:
        kwargs['language'] = 'fortran'
        super().__init__(name, environment, kwargs)

        cid = self.get_compiler().get_id()
        if cid == 'gcc':
            # Fallback to single image
            self.compile_args = ['-fcoarray=single']
            self.version = 'single image (fallback)'
            self.is_found = True
        elif cid == 'intel':
            # Coarrays are built into Intel compilers, no external library needed
            self.is_found = True
            self.link_args = ['-coarray=shared']
            self.compile_args = self.link_args
        elif cid == 'intel-cl':
            # Coarrays are built into Intel compilers, no external library needed
            self.is_found = True
            self.compile_args = ['/Qcoarray:shared']
        elif cid == 'nagfor':
            # NAG doesn't require any special arguments for Coarray
            self.is_found = True
