# SPDX-License-Identifier: Apache-2.0
# Copyright 2022 The Meson development team

from __future__ import annotations
import typing as T

from .backends import Backend
from .. import mlog
from ..mesonlib import MesonBugException


class NoneBackend(Backend):

    name = 'none'

    def generate(self, capture: bool = False, vslite_ctx: T.Optional[T.Dict] = None) -> None:
        # Check for (currently) unexpected capture arg use cases -
        if capture:
            raise MesonBugException('We do not expect the none backend to generate with \'capture = True\'')
        if vslite_ctx:
            raise MesonBugException('We do not expect the none backend to be given a valid \'vslite_ctx\'')

        if self.build.get_targets():
            raise MesonBugException('None backend cannot generate target rules, but should have failed earlier.')
        mlog.log('Generating simple install-only backend')
        self.serialize_tests()
        self.create_install_data_files()
