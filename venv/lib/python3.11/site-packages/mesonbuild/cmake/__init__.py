# SPDX-License-Identifier: Apache-2.0
# Copyright 2019 The Meson development team

# This class contains the basic functionality needed to run any interpreter
# or an interpreter-based tool.

__all__ = [
    'CMakeExecutor',
    'CMakeExecScope',
    'CMakeException',
    'CMakeInterpreter',
    'CMakeTarget',
    'CMakeToolchain',
    'CMakeTraceParser',
    'TargetOptions',
    'language_map',
    'cmake_defines_to_args',
    'check_cmake_args',
    'cmake_is_debug',
    'resolve_cmake_trace_targets',
]

from .common import CMakeException, TargetOptions, cmake_defines_to_args, language_map, check_cmake_args, cmake_is_debug
from .executor import CMakeExecutor
from .interpreter import CMakeInterpreter
from .toolchain import CMakeToolchain, CMakeExecScope
from .traceparser import CMakeTarget, CMakeTraceParser
from .tracetargets import resolve_cmake_trace_targets
