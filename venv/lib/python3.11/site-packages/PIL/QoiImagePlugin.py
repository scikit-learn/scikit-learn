#
# The Python Imaging Library.
#
# QOI support for PIL
#
# See the README file for information on usage and redistribution.
#
from __future__ import annotations

import os

from . import Image, ImageFile
from ._binary import i32be as i32


def _accept(prefix: bytes) -> bool:
    return prefix[:4] == b"qoif"


class QoiImageFile(ImageFile.ImageFile):
    format = "QOI"
    format_description = "Quite OK Image"

    def _open(self) -> None:
        if not _accept(self.fp.read(4)):
            msg = "not a QOI file"
            raise SyntaxError(msg)

        self._size = i32(self.fp.read(4)), i32(self.fp.read(4))

        channels = self.fp.read(1)[0]
        self._mode = "RGB" if channels == 3 else "RGBA"

        self.fp.seek(1, os.SEEK_CUR)  # colorspace
        self.tile = [ImageFile._Tile("qoi", (0, 0) + self._size, self.fp.tell())]


class QoiDecoder(ImageFile.PyDecoder):
    _pulls_fd = True
    _previous_pixel: bytes | bytearray | None = None
    _previously_seen_pixels: dict[int, bytes | bytearray] = {}

    def _add_to_previous_pixels(self, value: bytes | bytearray) -> None:
        self._previous_pixel = value

        r, g, b, a = value
        hash_value = (r * 3 + g * 5 + b * 7 + a * 11) % 64
        self._previously_seen_pixels[hash_value] = value

    def decode(self, buffer: bytes | Image.SupportsArrayInterface) -> tuple[int, int]:
        assert self.fd is not None

        self._previously_seen_pixels = {}
        self._add_to_previous_pixels(bytearray((0, 0, 0, 255)))

        data = bytearray()
        bands = Image.getmodebands(self.mode)
        dest_length = self.state.xsize * self.state.ysize * bands
        while len(data) < dest_length:
            byte = self.fd.read(1)[0]
            value: bytes | bytearray
            if byte == 0b11111110 and self._previous_pixel:  # QOI_OP_RGB
                value = bytearray(self.fd.read(3)) + self._previous_pixel[3:]
            elif byte == 0b11111111:  # QOI_OP_RGBA
                value = self.fd.read(4)
            else:
                op = byte >> 6
                if op == 0:  # QOI_OP_INDEX
                    op_index = byte & 0b00111111
                    value = self._previously_seen_pixels.get(
                        op_index, bytearray((0, 0, 0, 0))
                    )
                elif op == 1 and self._previous_pixel:  # QOI_OP_DIFF
                    value = bytearray(
                        (
                            (self._previous_pixel[0] + ((byte & 0b00110000) >> 4) - 2)
                            % 256,
                            (self._previous_pixel[1] + ((byte & 0b00001100) >> 2) - 2)
                            % 256,
                            (self._previous_pixel[2] + (byte & 0b00000011) - 2) % 256,
                            self._previous_pixel[3],
                        )
                    )
                elif op == 2 and self._previous_pixel:  # QOI_OP_LUMA
                    second_byte = self.fd.read(1)[0]
                    diff_green = (byte & 0b00111111) - 32
                    diff_red = ((second_byte & 0b11110000) >> 4) - 8
                    diff_blue = (second_byte & 0b00001111) - 8

                    value = bytearray(
                        tuple(
                            (self._previous_pixel[i] + diff_green + diff) % 256
                            for i, diff in enumerate((diff_red, 0, diff_blue))
                        )
                    )
                    value += self._previous_pixel[3:]
                elif op == 3 and self._previous_pixel:  # QOI_OP_RUN
                    run_length = (byte & 0b00111111) + 1
                    value = self._previous_pixel
                    if bands == 3:
                        value = value[:3]
                    data += value * run_length
                    continue
            self._add_to_previous_pixels(value)

            if bands == 3:
                value = value[:3]
            data += value
        self.set_as_raw(data)
        return -1, 0


Image.register_open(QoiImageFile.format, QoiImageFile, _accept)
Image.register_decoder("qoi", QoiDecoder)
Image.register_extension(QoiImageFile.format, ".qoi")
