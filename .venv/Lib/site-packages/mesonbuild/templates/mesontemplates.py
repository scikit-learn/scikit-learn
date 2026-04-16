# SPDX-License-Identifier: Apache-2.0
# Copyright 2019 The Meson development team
# Copyright © 2023-2025 Intel Corporation

from __future__ import annotations

import typing as T

from .samplefactory import sample_generator

if T.TYPE_CHECKING:
    from ..minit import Arguments


def create_meson_build(options: Arguments) -> None:
    proj = sample_generator(options)
    if options.type == 'executable':
        proj.create_executable()
    else:
        proj.create_library()
