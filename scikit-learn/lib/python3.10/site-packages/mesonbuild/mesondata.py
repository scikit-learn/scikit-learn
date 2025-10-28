# SPDX-License-Identifier: Apache-2.0
# Copyright 2021 The Meson development team

from __future__ import annotations


import importlib.resources
from pathlib import PurePosixPath, Path
import typing as T

if T.TYPE_CHECKING:
    from .environment import Environment

class DataFile:
    def __init__(self, path: str) -> None:
        self.path = PurePosixPath(path)

    def write_once(self, path: Path) -> None:
        if not path.exists():
            data = importlib.resources.read_text( # [ignore encoding] it's on the next lines, Mr. Lint
                    ('mesonbuild' / self.path.parent).as_posix().replace('/', '.'),
                    self.path.name,
                    encoding='utf-8')
            path.write_text(data, encoding='utf-8')

    def write_to_private(self, env: 'Environment') -> Path:
        try:
            resource = importlib.resources.files('mesonbuild') / self.path
            if isinstance(resource, Path):
                return resource
        except AttributeError:
            # fall through to python 3.7 compatible code
            pass

        out_file = Path(env.scratch_dir) / 'data' / self.path.name
        out_file.parent.mkdir(exist_ok=True)
        self.write_once(out_file)
        return out_file
