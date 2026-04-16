# SPDX-License-Identifier: Apache-2.0
# Copyright 2019 The Meson development team
# Copyright © 2023-2025 Intel Corporation

from __future__ import annotations

import typing as T

from mesonbuild.templates.valatemplates import ValaProject
from mesonbuild.templates.fortrantemplates import FortranProject
from mesonbuild.templates.objcpptemplates import ObjCppProject
from mesonbuild.templates.dlangtemplates import DlangProject
from mesonbuild.templates.rusttemplates import RustProject
from mesonbuild.templates.javatemplates import JavaProject
from mesonbuild.templates.cudatemplates import CudaProject
from mesonbuild.templates.objctemplates import ObjCProject
from mesonbuild.templates.cpptemplates import CppProject
from mesonbuild.templates.cstemplates import CSharpProject
from mesonbuild.templates.ctemplates import CProject

if T.TYPE_CHECKING:
    from ..minit import Arguments
    from .sampleimpl import ClassImpl, FileHeaderImpl, FileImpl, SampleImpl


_IMPL: T.Mapping[str, T.Union[T.Type[ClassImpl], T.Type[FileHeaderImpl], T.Type[FileImpl]]] = {
    'c': CProject,
    'cpp': CppProject,
    'cs': CSharpProject,
    'cuda': CudaProject,
    'objc': ObjCProject,
    'objcpp': ObjCppProject,
    'java': JavaProject,
    'd': DlangProject,
    'rust': RustProject,
    'fortran': FortranProject,
    'vala': ValaProject,
}


def sample_generator(options: Arguments) -> SampleImpl:
    return _IMPL[options.language](options)
