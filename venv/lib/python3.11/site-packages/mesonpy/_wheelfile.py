# SPDX-FileCopyrightText: 2022 The meson-python developers
#
# SPDX-License-Identifier: MIT

from __future__ import annotations

import base64
import csv
import hashlib
import io
import os
import re
import stat
import time
import typing
import zipfile


if typing.TYPE_CHECKING:  # pragma: no cover
    from types import TracebackType
    from typing import List, Optional, Tuple, Type, Union

    from mesonpy._compat import Path


MIN_TIMESTAMP = 315532800  # 1980-01-01 00:00:00 UTC
WHEEL_FILENAME_REGEX = re.compile(r'^(?P<name>[^-]+)-(?P<version>[^-]+)(:?-(?P<build>[^-]+))?-(?P<tag>[^-]+-[^-]+-[^-]+).whl$')


def _b64encode(data: bytes) -> bytes:
    return base64.urlsafe_b64encode(data).rstrip(b'=')


class WheelFile:
    """Implement the wheel package binary distribution format.

    https://packaging.python.org/en/latest/specifications/binary-distribution-format/
    """
    def __new__(cls, filename: Path, mode: str = 'r', compression: int = zipfile.ZIP_DEFLATED) -> 'WheelFile':
        if mode == 'w':
            return super().__new__(WheelFileWriter)
        raise NotImplementedError

    @staticmethod
    def timestamp(mtime: Optional[float] = None) -> Tuple[int, int, int, int, int, int]:
        timestamp = int(os.environ.get('SOURCE_DATE_EPOCH', mtime or time.time()))
        # The ZIP file format does not support timestamps before 1980.
        timestamp = max(timestamp, MIN_TIMESTAMP)
        return time.gmtime(timestamp)[0:6]

    @staticmethod
    def hash(data: bytes) -> str:
        return 'sha256=' + _b64encode(hashlib.sha256(data).digest()).decode('ascii')

    def writestr(self, zinfo_or_arcname: Union[str, zipfile.ZipInfo], data: bytes) -> None:
        raise NotImplementedError

    def write(self, filename: Path, arcname: Optional[str] = None) -> None:
        raise NotImplementedError

    def close(self) -> None:
        raise NotImplementedError

    def __enter__(self) -> WheelFile:
        return self

    def __exit__(self, exc_type: Type[BaseException], exc_val: BaseException, exc_tb: TracebackType) -> None:
        self.close()


class WheelFileWriter(WheelFile):
    def __init__(self, filepath: Path, mode: str, compression: int = zipfile.ZIP_DEFLATED):
        filename = os.path.basename(filepath)
        match = WHEEL_FILENAME_REGEX.match(filename)
        if not match:
            raise ValueError(f'invalid wheel filename: {filename!r}')
        self.name = match.group('name')
        self.version = match.group('version')
        self.entries: List[Tuple[str, str, int]] = []
        self.archive = zipfile.ZipFile(filepath, mode='w', compression=compression, allowZip64=True)

    def writestr(self, zinfo_or_arcname: Union[str, zipfile.ZipInfo], data: bytes) -> None:
        if isinstance(data, str):
            data = data.encode('utf-8')
        if isinstance(zinfo_or_arcname, zipfile.ZipInfo):
            zinfo = zinfo_or_arcname
        else:
            zinfo = zipfile.ZipInfo(zinfo_or_arcname, date_time=self.timestamp())
            zinfo.external_attr = 0o664 << 16
        self.archive.writestr(
            zinfo, data,
            compress_type=self.archive.compression,
            compresslevel=self.archive.compresslevel)
        self.entries.append((zinfo.filename, self.hash(data), len(data)))

    def write(self, filename: Path, arcname: Optional[str] = None) -> None:
        with open(filename, 'rb') as f:
            st = os.fstat(f.fileno())
            data = f.read()
        zinfo = zipfile.ZipInfo(arcname or str(filename), date_time=self.timestamp(st.st_mtime))
        zinfo.external_attr = (stat.S_IMODE(st.st_mode) | stat.S_IFMT(st.st_mode)) << 16
        self.writestr(zinfo, data)

    def close(self) -> None:
        record = f'{self.name}-{self.version}.dist-info/RECORD'
        data = io.StringIO()
        writer = csv.writer(data, delimiter=',', quotechar='"', lineterminator='\n')
        writer.writerows(self.entries)
        writer.writerow((record, '', ''))
        zi = zipfile.ZipInfo(record, date_time=self.timestamp())
        zi.external_attr = 0o664 << 16
        self.archive.writestr(
            zi, data.getvalue(),
            compress_type=self.archive.compression,
            compresslevel=self.archive.compresslevel)
        self.archive.close()
