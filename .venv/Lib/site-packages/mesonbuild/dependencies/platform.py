# SPDX-License-Identifier: Apache-2.0
# Copyright 2013-2017 The Meson development team

# This file contains the detection logic for external dependencies that are
# platform-specific (generally speaking).
from __future__ import annotations

from .base import DependencyTypeName, ExternalDependency, DependencyException
from .detect import packages
from ..mesonlib import MesonException
import typing as T

if T.TYPE_CHECKING:
    from ..environment import Environment
    from .base import DependencyObjectKWs

class AppleFrameworks(ExternalDependency):

    type_name = DependencyTypeName('appleframeworks')

    def __init__(self, name: str, env: 'Environment', kwargs: DependencyObjectKWs) -> None:
        super().__init__(name, env, kwargs)
        modules = kwargs.get('modules', [])
        if not modules:
            raise DependencyException("AppleFrameworks dependency requires at least one module.")
        self.frameworks = modules
        if not self.clib_compiler:
            raise DependencyException('No C-like compilers are available, cannot find the framework')
        self.is_found = True
        for f in self.frameworks:
            try:
                args = self.clib_compiler.find_framework(f, [])
            except MesonException as e:
                if 'non-clang' in str(e):
                    self.is_found = False
                    self.link_args = []
                    self.compile_args = []
                    return
                raise

            if args is not None:
                # No compile args are needed for system frameworks
                self.link_args += args
            else:
                self.is_found = False

    def log_info(self) -> str:
        return ', '.join(self.frameworks)

packages['appleframeworks'] = AppleFrameworks
