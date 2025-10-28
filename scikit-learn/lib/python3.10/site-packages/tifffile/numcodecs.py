# tifffile/numcodecs.py

# Copyright (c) 2021-2025, Christoph Gohlke
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""TIFF codec for the Numcodecs package."""

from __future__ import annotations

__all__ = ['register_codec', 'Tiff']

from io import BytesIO
from typing import TYPE_CHECKING

from numcodecs import registry
from numcodecs.abc import Codec

from .tifffile import TiffFile, TiffWriter

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence
    from typing import Any

    from .tifffile import (
        COMPRESSION,
        EXTRASAMPLE,
        PHOTOMETRIC,
        PLANARCONFIG,
        PREDICTOR,
        ByteOrder,
        TagTuple,
    )


class Tiff(Codec):  # type: ignore[misc]
    """TIFF codec for Numcodecs."""

    codec_id = 'tifffile'

    def __init__(
        self,
        # TiffFile.asarray
        key: int | slice | Iterable[int] | None = None,
        series: int | None = None,
        level: int | None = None,
        # TiffWriter
        bigtiff: bool = False,
        byteorder: ByteOrder | None = None,
        imagej: bool = False,
        ome: bool | None = None,
        # TiffWriter.write
        photometric: PHOTOMETRIC | int | str | None = None,
        planarconfig: PLANARCONFIG | int | str | None = None,
        extrasamples: Sequence[EXTRASAMPLE | int | str] | None = None,
        volumetric: bool = False,
        tile: Sequence[int] | None = None,
        rowsperstrip: int | None = None,
        compression: COMPRESSION | int | str | None = None,
        compressionargs: dict[str, Any] | None = None,
        predictor: PREDICTOR | int | str | bool | None = None,
        subsampling: tuple[int, int] | None = None,
        metadata: dict[str, Any] | None = {},
        extratags: Sequence[TagTuple] | None = None,
        truncate: bool = False,
        maxworkers: int | None = None,
    ) -> None:
        self.key = key
        self.series = series
        self.level = level
        self.bigtiff = bigtiff
        self.byteorder = byteorder
        self.imagej = imagej
        self.ome = ome
        self.photometric = photometric
        self.planarconfig = planarconfig
        self.extrasamples = extrasamples
        self.volumetric = volumetric
        self.tile = tile
        self.rowsperstrip = rowsperstrip
        self.compression = compression
        self.compressionargs = compressionargs
        self.predictor = predictor
        self.subsampling = subsampling
        self.metadata = metadata
        self.extratags = extratags
        self.truncate = truncate
        self.maxworkers = maxworkers

    def encode(self, buf: Any) -> bytes:
        """Return TIFF file as bytes."""
        with BytesIO() as fh:
            with TiffWriter(
                fh,
                bigtiff=self.bigtiff,
                byteorder=self.byteorder,
                imagej=self.imagej,
                ome=self.ome,
            ) as tif:
                tif.write(
                    buf,
                    photometric=self.photometric,
                    planarconfig=self.planarconfig,
                    extrasamples=self.extrasamples,
                    volumetric=self.volumetric,
                    tile=self.tile,
                    rowsperstrip=self.rowsperstrip,
                    compression=self.compression,
                    compressionargs=self.compressionargs,
                    predictor=self.predictor,
                    subsampling=self.subsampling,
                    metadata=self.metadata,
                    extratags=self.extratags,
                    truncate=self.truncate,
                    maxworkers=self.maxworkers,
                )
            result = fh.getvalue()
        return result

    def decode(self, buf: Any, out: Any = None) -> Any:
        """Return decoded image as NumPy array."""
        with BytesIO(buf) as fh:
            with TiffFile(fh) as tif:
                result = tif.asarray(
                    key=self.key,
                    series=self.series,
                    level=self.level,
                    maxworkers=self.maxworkers,
                    out=out,
                )
        return result


def register_codec(cls: Codec = Tiff, codec_id: str | None = None) -> None:
    """Register :py:class:`Tiff` codec with Numcodecs."""
    registry.register_codec(cls, codec_id=codec_id)
