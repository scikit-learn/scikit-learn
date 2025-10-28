# -*- coding: utf-8 -*-
# imageio is distributed under the terms of the (new) BSD License.

""" Read SPE files.

This plugin supports reading files saved in the Princeton Instruments
SPE file format.

Parameters
----------
check_filesize : bool
    The number of frames in the file is stored in the file header. However,
    this number may be wrong for certain software. If this is `True`
    (default), derive the number of frames also from the file size and
    raise a warning if the two values do not match.
char_encoding : str
    Deprecated. Exists for backwards compatibility; use ``char_encoding`` of
    ``metadata`` instead.
sdt_meta : bool
    Deprecated. Exists for backwards compatibility; use ``sdt_control`` of
    ``metadata`` instead.

Methods
-------
.. note::
    Check the respective function for a list of supported kwargs and detailed
    documentation.

.. autosummary::
    :toctree:

    SpePlugin.read
    SpePlugin.iter
    SpePlugin.properties
    SpePlugin.metadata

"""

from datetime import datetime
import logging
import os
from typing import (
    Any,
    Callable,
    Dict,
    Iterator,
    List,
    Mapping,
    Optional,
    Sequence,
    Tuple,
    Union,
)
import warnings

import numpy as np

from ..core.request import Request, IOMode, InitializationError
from ..core.v3_plugin_api import PluginV3, ImageProperties


logger = logging.getLogger(__name__)


class Spec:
    """SPE file specification data

    Tuples of (offset, datatype, count), where offset is the offset in the SPE
    file and datatype is the datatype as used in `numpy.fromfile`()

    `data_start` is the offset of actual image data.

    `dtypes` translates SPE datatypes (0...4) to numpy ones, e. g. dtypes[0]
    is dtype("<f") (which is np.float32).

    `controllers` maps the `type` metadata to a human readable name

    `readout_modes` maps the `readoutMode` metadata to something human readable
    although this may not be accurate since there is next to no documentation
    to be found.
    """

    basic = {
        "datatype": (108, "<h"),  # dtypes
        "xdim": (42, "<H"),
        "ydim": (656, "<H"),
        "xml_footer_offset": (678, "<Q"),
        "NumFrames": (1446, "<i"),
        "file_header_ver": (1992, "<f"),
    }

    metadata = {
        # ROI information
        "NumROI": (1510, "<h"),
        "ROIs": (
            1512,
            np.dtype(
                [
                    ("startx", "<H"),
                    ("endx", "<H"),
                    ("groupx", "<H"),
                    ("starty", "<H"),
                    ("endy", "<H"),
                    ("groupy", "<H"),
                ]
            ),
            10,
        ),
        # chip-related sizes
        "xDimDet": (6, "<H"),
        "yDimDet": (18, "<H"),
        "VChipXdim": (14, "<h"),
        "VChipYdim": (16, "<h"),
        # other stuff
        "controller_version": (0, "<h"),
        "logic_output": (2, "<h"),
        "amp_high_cap_low_noise": (4, "<H"),  # enum?
        "mode": (8, "<h"),  # enum?
        "exposure_sec": (10, "<f"),
        "date": (20, "<10S"),
        "detector_temp": (36, "<f"),
        "detector_type": (40, "<h"),
        "st_diode": (44, "<h"),
        "delay_time": (46, "<f"),
        # shutter_control: normal, disabled open, disabled closed
        # But which one is which?
        "shutter_control": (50, "<H"),
        "absorb_live": (52, "<h"),
        "absorb_mode": (54, "<H"),
        "can_do_virtual_chip": (56, "<h"),
        "threshold_min_live": (58, "<h"),
        "threshold_min_val": (60, "<f"),
        "threshold_max_live": (64, "<h"),
        "threshold_max_val": (66, "<f"),
        "time_local": (172, "<7S"),
        "time_utc": (179, "<7S"),
        "adc_offset": (188, "<H"),
        "adc_rate": (190, "<H"),
        "adc_type": (192, "<H"),
        "adc_resolution": (194, "<H"),
        "adc_bit_adjust": (196, "<H"),
        "gain": (198, "<H"),
        "comments": (200, "<80S", 5),
        "geometric": (600, "<H"),  # flags
        "sw_version": (688, "<16S"),
        "spare_4": (742, "<436S"),
        "XPrePixels": (98, "<h"),
        "XPostPixels": (100, "<h"),
        "YPrePixels": (102, "<h"),
        "YPostPixels": (104, "<h"),
        "readout_time": (672, "<f"),
        "xml_footer_offset": (678, "<Q"),
        "type": (704, "<h"),  # controllers
        "clockspeed_us": (1428, "<f"),
        "readout_mode": (1480, "<H"),  # readout_modes
        "window_size": (1482, "<H"),
        "file_header_ver": (1992, "<f"),
    }

    data_start = 4100

    dtypes = {
        0: np.dtype(np.float32),
        1: np.dtype(np.int32),
        2: np.dtype(np.int16),
        3: np.dtype(np.uint16),
        8: np.dtype(np.uint32),
    }

    controllers = [
        "new120 (Type II)",
        "old120 (Type I)",
        "ST130",
        "ST121",
        "ST138",
        "DC131 (PentaMax)",
        "ST133 (MicroMax/Roper)",
        "ST135 (GPIB)",
        "VTCCD",
        "ST116 (GPIB)",
        "OMA3 (GPIB)",
        "OMA4",
    ]

    # This was gathered from random places on the internet and own experiments
    # with the camera. May not be accurate.
    readout_modes = ["full frame", "frame transfer", "kinetics"]

    # Do not decode the following metadata keys into strings, but leave them
    # as byte arrays
    no_decode = ["spare_4"]


class SDTControlSpec:
    """Extract metadata written by the SDT-control software

    Some of it is encoded in the comment strings
    (see :py:meth:`parse_comments`). Also, date and time are encoded in a
    peculiar way (see :py:meth:`get_datetime`). Use :py:meth:`extract_metadata`
    to update the metadata dict.
    """

    months = {
        # Convert SDT-control month strings to month numbers
        "Jän": 1,
        "Jan": 1,
        "Feb": 2,
        "Mär": 3,
        "Mar": 3,
        "Apr": 4,
        "Mai": 5,
        "May": 5,
        "Jun": 6,
        "Jul": 7,
        "Aug": 8,
        "Sep": 9,
        "Okt": 10,
        "Oct": 10,
        "Nov": 11,
        "Dez": 12,
        "Dec": 12,
    }

    sequence_types = {
        # TODO: complete
        "SEQU": "standard",
        "SETO": "TOCCSL",
        "KINE": "kinetics",
        "SEAR": "arbitrary",
    }

    class CommentDesc:
        """Describe how to extract a metadata entry from a comment string"""

        n: int
        """Which of the 5 SPE comment fields to use."""
        slice: slice
        """Which characters from the `n`-th comment to use."""
        cvt: Callable[[str], Any]
        """How to convert characters to something useful."""
        scale: Union[None, float]
        """Optional scaling factor for numbers"""

        def __init__(
            self,
            n: int,
            slice: slice,
            cvt: Callable[[str], Any] = str,
            scale: Optional[float] = None,
        ):
            self.n = n
            self.slice = slice
            self.cvt = cvt
            self.scale = scale

    comment_fields = {
        (5, 0): {
            "sdt_major_version": CommentDesc(4, slice(66, 68), int),
            "sdt_minor_version": CommentDesc(4, slice(68, 70), int),
            "sdt_controller_name": CommentDesc(4, slice(0, 6), str),
            "exposure_time": CommentDesc(1, slice(64, 73), float, 10**-6),
            "color_code": CommentDesc(4, slice(10, 14), str),
            "detection_channels": CommentDesc(4, slice(15, 16), int),
            "background_subtraction": CommentDesc(4, 14, lambda x: x == "B"),
            "em_active": CommentDesc(4, 32, lambda x: x == "E"),
            "em_gain": CommentDesc(4, slice(28, 32), int),
            "modulation_active": CommentDesc(4, 33, lambda x: x == "A"),
            "pixel_size": CommentDesc(4, slice(25, 28), float, 0.1),
            "sequence_type": CommentDesc(
                4, slice(6, 10), lambda x: __class__.sequence_types[x]
            ),
            "grid": CommentDesc(4, slice(16, 25), float, 10**-6),
            "n_macro": CommentDesc(1, slice(0, 4), int),
            "delay_macro": CommentDesc(1, slice(10, 19), float, 10**-3),
            "n_mini": CommentDesc(1, slice(4, 7), int),
            "delay_mini": CommentDesc(1, slice(19, 28), float, 10**-6),
            "n_micro": CommentDesc(1, slice(7, 10), int),
            "delay_micro": CommentDesc(1, slice(28, 37), float, 10**-6),
            "n_subpics": CommentDesc(1, slice(7, 10), int),
            "delay_shutter": CommentDesc(1, slice(73, 79), float, 10**-6),
            "delay_prebleach": CommentDesc(1, slice(37, 46), float, 10**-6),
            "bleach_time": CommentDesc(1, slice(46, 55), float, 10**-6),
            "recovery_time": CommentDesc(1, slice(55, 64), float, 10**-6),
        },
        (5, 1): {
            "bleach_piezo_active": CommentDesc(4, slice(34, 35), lambda x: x == "z")
        },
    }

    @staticmethod
    def get_comment_version(comments: Sequence[str]) -> Tuple[int, int]:
        """Get the version of SDT-control metadata encoded in the comments

        Parameters
        ----------
        comments
            List of SPE file comments, typically ``metadata["comments"]``.

        Returns
        -------
        Major and minor version. ``-1, -1`` if detection failed.
        """
        if comments[4][70:76] != "COMVER":
            return -1, -1
        try:
            return int(comments[4][76:78]), int(comments[4][78:80])
        except ValueError:
            return -1, -1

    @staticmethod
    def parse_comments(
        comments: Sequence[str], version: Tuple[int, int]
    ) -> Dict[str, Any]:
        """Extract SDT-control metadata from comments

        Parameters
        ----------
        comments
            List of SPE file comments, typically ``metadata["comments"]``.
        version
            Major and minor version of SDT-control metadata format

        Returns
        -------
        Dict of metadata
        """
        sdt_md = {}
        for minor in range(version[1] + 1):
            # Metadata with same major version is backwards compatible.
            # Fields are specified incrementally in `comment_fields`.
            # E.g. if the file has version 5.01, `comment_fields[5, 0]` and
            # `comment_fields[5, 1]` need to be decoded.
            try:
                cmt = __class__.comment_fields[version[0], minor]
            except KeyError:
                continue
            for name, spec in cmt.items():
                try:
                    v = spec.cvt(comments[spec.n][spec.slice])
                    if spec.scale is not None:
                        v *= spec.scale
                    sdt_md[name] = v
                except Exception as e:
                    warnings.warn(
                        f"Failed to decode SDT-control metadata field `{name}`: {e}"
                    )
                    sdt_md[name] = None
        if version not in __class__.comment_fields:
            supported_ver = ", ".join(
                map(lambda x: f"{x[0]}.{x[1]:02}", __class__.comment_fields)
            )
            warnings.warn(
                f"Unsupported SDT-control metadata version {version[0]}.{version[1]:02}. "
                f"Only versions {supported_ver} are supported. "
                "Some or all SDT-control metadata may be missing."
            )
        comment = comments[0] + comments[2]
        sdt_md["comment"] = comment.strip()
        return sdt_md

    @staticmethod
    def get_datetime(date: str, time: str) -> Union[datetime, None]:
        """Turn date and time saved by SDT-control into proper datetime object

        Parameters
        ----------
        date
            SPE file date, typically ``metadata["date"]``.
        time
            SPE file date, typically ``metadata["time_local"]``.

        Returns
        -------
        File's datetime if parsing was succsessful, else None.
        """
        try:
            month = __class__.months[date[2:5]]
            return datetime(
                int(date[5:9]),
                month,
                int(date[0:2]),
                int(time[0:2]),
                int(time[2:4]),
                int(time[4:6]),
            )
        except Exception as e:
            logger.info(f"Failed to decode date from SDT-control metadata: {e}.")

    @staticmethod
    def extract_metadata(meta: Mapping, char_encoding: str = "latin1"):
        """Extract SDT-control metadata from SPE metadata

        SDT-control stores some metadata in comments and other fields.
        Extract them and remove unused entries.

        Parameters
        ----------
        meta
            SPE file metadata. Modified in place.
        char_encoding
            Character encoding used to decode strings in the metadata.
        """
        comver = __class__.get_comment_version(meta["comments"])
        if any(c < 0 for c in comver):
            # This file most likely was not created by SDT-control
            logger.debug("SDT-control comments not found.")
            return

        sdt_meta = __class__.parse_comments(meta["comments"], comver)
        meta.pop("comments")
        meta.update(sdt_meta)

        # Get date and time in a usable format
        dt = __class__.get_datetime(meta["date"], meta["time_local"])
        if dt:
            meta["datetime"] = dt
            meta.pop("date")
            meta.pop("time_local")

        sp4 = meta["spare_4"]
        try:
            meta["modulation_script"] = sp4.decode(char_encoding)
            meta.pop("spare_4")
        except UnicodeDecodeError:
            warnings.warn(
                "Failed to decode SDT-control laser "
                "modulation script. Bad char_encoding?"
            )

        # Get rid of unused data
        meta.pop("time_utc")
        meta.pop("exposure_sec")


class SpePlugin(PluginV3):
    def __init__(
        self,
        request: Request,
        check_filesize: bool = True,
        char_encoding: Optional[str] = None,
        sdt_meta: Optional[bool] = None,
    ) -> None:
        """Instantiate a new SPE file plugin object

        Parameters
        ----------
        request : Request
            A request object representing the resource to be operated on.
        check_filesize : bool
            If True, compute the number of frames from the filesize, compare it
            to the frame count in the file header, and raise a warning if the
            counts don't match. (Certain software may create files with
        char_encoding : str
            Deprecated. Exists for backwards compatibility; use ``char_encoding`` of
            ``metadata`` instead.
        sdt_meta : bool
            Deprecated. Exists for backwards compatibility; use ``sdt_control`` of
            ``metadata`` instead.

        """

        super().__init__(request)
        if request.mode.io_mode == IOMode.write:
            raise InitializationError("cannot write SPE files")

        if char_encoding is not None:
            warnings.warn(
                "Passing `char_encoding` to the constructor is deprecated. "
                "Use `char_encoding` parameter of the `metadata()` method "
                "instead.",
                DeprecationWarning,
            )
        self._char_encoding = char_encoding
        if sdt_meta is not None:
            warnings.warn(
                "Passing `sdt_meta` to the constructor is deprecated. "
                "Use `sdt_control` parameter of the `metadata()` method "
                "instead.",
                DeprecationWarning,
            )
        self._sdt_meta = sdt_meta

        self._file = self.request.get_file()

        try:
            # Spec.basic contains no string, no need to worry about character
            # encoding.
            info = self._parse_header(Spec.basic, "latin1")
            self._file_header_ver = info["file_header_ver"]
            self._dtype = Spec.dtypes[info["datatype"]]
            self._shape = (info["ydim"], info["xdim"])
            self._len = info["NumFrames"]

            if check_filesize:
                # Some software writes incorrect `NumFrames` metadata.
                # To determine the number of frames, check the size of the data
                # segment -- until the end of the file for SPE<3, until the
                # xml footer for SPE>=3.
                if info["file_header_ver"] >= 3:
                    data_end = info["xml_footer_offset"]
                else:
                    self._file.seek(0, os.SEEK_END)
                    data_end = self._file.tell()
                line = data_end - Spec.data_start
                line //= self._shape[0] * self._shape[1] * self._dtype.itemsize
                if line != self._len:
                    warnings.warn(
                        f"The file header of {self.request.filename} claims there are "
                        f"{self._len} frames, but there are actually {line} frames."
                    )
                    self._len = min(line, self._len)
            self._file.seek(Spec.data_start)
        except Exception:
            raise InitializationError("SPE plugin cannot read the provided file.")

    def read(self, *, index: int = ...) -> np.ndarray:
        """Read a frame or all frames from the file

        Parameters
        ----------
        index : int
            Select the index-th frame from the file. If index is `...`,
            select all frames and stack them along a new axis.

        Returns
        -------
        A Numpy array of pixel values.

        """

        if index is Ellipsis:
            read_offset = Spec.data_start
            count = self._shape[0] * self._shape[1] * self._len
            out_shape = (self._len, *self._shape)
        elif index < 0:
            raise IndexError(f"Index `{index}` is smaller than 0.")
        elif index >= self._len:
            raise IndexError(
                f"Index `{index}` exceeds the number of frames stored in this file (`{self._len}`)."
            )
        else:
            read_offset = (
                Spec.data_start
                + index * self._shape[0] * self._shape[1] * self._dtype.itemsize
            )
            count = self._shape[0] * self._shape[1]
            out_shape = self._shape

        self._file.seek(read_offset)
        data = np.fromfile(self._file, dtype=self._dtype, count=count)
        return data.reshape(out_shape)

    def iter(self) -> Iterator[np.ndarray]:
        """Iterate over the frames in the file

        Yields
        ------
        A Numpy array of pixel values.
        """

        return (self.read(index=i) for i in range(self._len))

    def metadata(
        self,
        index: int = ...,
        exclude_applied: bool = True,
        char_encoding: str = "latin1",
        sdt_control: bool = True,
    ) -> Dict[str, Any]:
        """SPE specific metadata.

        Parameters
        ----------
        index : int
            Ignored as SPE files only store global metadata.
        exclude_applied : bool
            Ignored. Exists for API compatibility.
        char_encoding : str
            The encoding to use when parsing strings.
        sdt_control : bool
            If `True`, decode special metadata written by the
            SDT-control software if present.

        Returns
        -------
        metadata : dict
            Key-value pairs of metadata.

        Notes
        -----
        SPE v3 stores metadata as XML, whereas SPE v2 uses a binary format.

        .. rubric:: Supported SPE v2 Metadata fields

        ROIs : list of dict
            Regions of interest used for recording images. Each dict has the
            "top_left" key containing x and y coordinates of the top left corner,
            the "bottom_right" key with x and y coordinates of the bottom right
            corner, and the "bin" key with number of binned pixels in x and y
            directions.
        comments : list of str
            The SPE format allows for 5 comment strings of 80 characters each.
        controller_version : int
            Hardware version
        logic_output : int
            Definition of output BNC
        amp_hi_cap_low_noise : int
            Amp switching mode
        mode : int
            Timing mode
        exp_sec : float
            Alternative exposure in seconds
        date : str
            Date string
        detector_temp : float
            Detector temperature
        detector_type : int
            CCD / diode array type
        st_diode : int
            Trigger diode
        delay_time : float
            Used with async mode
        shutter_control : int
            Normal, disabled open, or disabled closed
        absorb_live : bool
            on / off
        absorb_mode : int
            Reference strip or file
        can_do_virtual_chip : bool
            True or False whether chip can do virtual chip
        threshold_min_live : bool
            on / off
        threshold_min_val : float
            Threshold minimum value
        threshold_max_live : bool
            on / off
        threshold_max_val : float
            Threshold maximum value
        time_local : str
            Experiment local time
        time_utc : str
            Experiment UTC time
        adc_offset : int
            ADC offset
        adc_rate : int
            ADC rate
        adc_type : int
            ADC type
        adc_resolution : int
            ADC resolution
        adc_bit_adjust : int
            ADC bit adjust
        gain : int
            gain
        sw_version : str
            Version of software which created this file
        spare_4 : bytes
            Reserved space
        readout_time : float
            Experiment readout time
        type : str
            Controller type
        clockspeed_us : float
            Vertical clock speed in microseconds
        readout_mode : ["full frame", "frame transfer", "kinetics", ""]
            Readout mode. Empty string means that this was not set by the
            Software.
        window_size : int
            Window size for Kinetics mode
        file_header_ver : float
            File header version
        chip_size : [int, int]
            x and y dimensions of the camera chip
        virt_chip_size : [int, int]
            Virtual chip x and y dimensions
        pre_pixels : [int, int]
            Pre pixels in x and y dimensions
        post_pixels : [int, int],
            Post pixels in x and y dimensions
        geometric : list of {"rotate", "reverse", "flip"}
            Geometric operations
        sdt_major_version : int
            (only for files created by SDT-control)
            Major version of SDT-control software
        sdt_minor_version : int
            (only for files created by SDT-control)
            Minor version of SDT-control software
        sdt_controller_name : str
            (only for files created by SDT-control)
            Controller name
        exposure_time : float
            (only for files created by SDT-control)
            Exposure time in seconds
        color_code : str
            (only for files created by SDT-control)
            Color channels used
        detection_channels : int
            (only for files created by SDT-control)
            Number of channels
        background_subtraction : bool
            (only for files created by SDT-control)
            Whether background subtraction war turned on
        em_active : bool
            (only for files created by SDT-control)
            Whether EM was turned on
        em_gain : int
            (only for files created by SDT-control)
            EM gain
        modulation_active : bool
            (only for files created by SDT-control)
            Whether laser modulation (“attenuate”) was turned on
        pixel_size : float
            (only for files created by SDT-control)
            Camera pixel size
        sequence_type : str
            (only for files created by SDT-control)
            Type of sequnce (standard, TOCCSL, arbitrary, …)
        grid : float
            (only for files created by SDT-control)
            Sequence time unit (“grid size”) in seconds
        n_macro : int
            (only for files created by SDT-control)
            Number of macro loops
        delay_macro : float
            (only for files created by SDT-control)
            Time between macro loops in seconds
        n_mini : int
            (only for files created by SDT-control)
            Number of mini loops
        delay_mini : float
            (only for files created by SDT-control)
            Time between mini loops in seconds
        n_micro : int (only for files created by SDT-control)
            Number of micro loops
        delay_micro : float (only for files created by SDT-control)
            Time between micro loops in seconds
        n_subpics : int
            (only for files created by SDT-control)
            Number of sub-pictures
        delay_shutter : float
            (only for files created by SDT-control)
            Camera shutter delay in seconds
        delay_prebleach : float
            (only for files created by SDT-control)
            Pre-bleach delay in seconds
        bleach_time : float
            (only for files created by SDT-control)
            Bleaching time in seconds
        recovery_time : float
            (only for files created by SDT-control)
            Recovery time in seconds
        comment : str
            (only for files created by SDT-control)
            User-entered comment. This replaces the "comments" field.
        datetime : datetime.datetime
            (only for files created by SDT-control)
            Combines the "date" and "time_local" keys. The latter two plus
            "time_utc" are removed.
        modulation_script : str
            (only for files created by SDT-control)
            Laser modulation script. Replaces the "spare_4" key.
        bleach_piezo_active : bool
            (only for files created by SDT-control)
            Whether piezo for bleaching was enabled
        """

        if self._file_header_ver < 3:
            if self._char_encoding is not None:
                char_encoding = self._char_encoding
            if self._sdt_meta is not None:
                sdt_control = self._sdt_meta
            return self._metadata_pre_v3(char_encoding, sdt_control)
        return self._metadata_post_v3()

    def _metadata_pre_v3(self, char_encoding: str, sdt_control: bool) -> Dict[str, Any]:
        """Extract metadata from SPE v2 files

        Parameters
        ----------
        char_encoding
            String character encoding
        sdt_control
            If `True`, try to decode special metadata written by the
            SDT-control software.

        Returns
        -------
        dict mapping metadata names to values.

        """

        m = self._parse_header(Spec.metadata, char_encoding)

        nr = m.pop("NumROI", None)
        nr = 1 if nr < 1 else nr
        m["ROIs"] = roi_array_to_dict(m["ROIs"][:nr])

        # chip sizes
        m["chip_size"] = [m.pop(k, None) for k in ("xDimDet", "yDimDet")]
        m["virt_chip_size"] = [m.pop(k, None) for k in ("VChipXdim", "VChipYdim")]
        m["pre_pixels"] = [m.pop(k, None) for k in ("XPrePixels", "YPrePixels")]
        m["post_pixels"] = [m.pop(k, None) for k in ("XPostPixels", "YPostPixels")]

        # convert comments from numpy.str_ to str
        m["comments"] = [str(c) for c in m["comments"]]

        # geometric operations
        g = []
        f = m.pop("geometric", 0)
        if f & 1:
            g.append("rotate")
        if f & 2:
            g.append("reverse")
        if f & 4:
            g.append("flip")
        m["geometric"] = g

        # Make some additional information more human-readable
        t = m["type"]
        if 1 <= t <= len(Spec.controllers):
            m["type"] = Spec.controllers[t - 1]
        else:
            m["type"] = None
        r = m["readout_mode"]
        if 1 <= r <= len(Spec.readout_modes):
            m["readout_mode"] = Spec.readout_modes[r - 1]
        else:
            m["readout_mode"] = None

        # bools
        for k in (
            "absorb_live",
            "can_do_virtual_chip",
            "threshold_min_live",
            "threshold_max_live",
        ):
            m[k] = bool(m[k])

        # Extract SDT-control metadata if desired
        if sdt_control:
            SDTControlSpec.extract_metadata(m, char_encoding)

        return m

    def _metadata_post_v3(self) -> Dict[str, Any]:
        """Extract XML metadata from SPE v3 files

        Returns
        -------
        dict with key `"__xml"`, whose value is the XML metadata
        """

        info = self._parse_header(Spec.basic, "latin1")
        self._file.seek(info["xml_footer_offset"])
        xml = self._file.read()
        return {"__xml": xml}

    def properties(self, index: int = ...) -> ImageProperties:
        """Standardized ndimage metadata.

        Parameters
        ----------
        index : int
            If the index is an integer, select the index-th frame and return
            its properties. If index is an Ellipsis (...), return the
            properties of all frames in the file stacked along a new batch
            dimension.

        Returns
        -------
        properties : ImageProperties
            A dataclass filled with standardized image metadata.
        """

        if index is Ellipsis:
            return ImageProperties(
                shape=(self._len, *self._shape),
                dtype=self._dtype,
                n_images=self._len,
                is_batch=True,
            )
        return ImageProperties(shape=self._shape, dtype=self._dtype, is_batch=False)

    def _parse_header(
        self, spec: Mapping[str, Tuple], char_encoding: str
    ) -> Dict[str, Any]:
        """Get information from SPE file header

        Parameters
        ----------
        spec
            Maps header entry name to its location, data type description and
            optionally number of entries. See :py:attr:`Spec.basic` and
            :py:attr:`Spec.metadata`.
        char_encoding
            String character encoding

        Returns
        -------
        Dict mapping header entry name to its value
        """

        ret = {}
        # Decode each string from the numpy array read by np.fromfile
        decode = np.vectorize(lambda x: x.decode(char_encoding))

        for name, sp in spec.items():
            self._file.seek(sp[0])
            cnt = 1 if len(sp) < 3 else sp[2]
            v = np.fromfile(self._file, dtype=sp[1], count=cnt)
            if v.dtype.kind == "S" and name not in Spec.no_decode:
                # Silently ignore string decoding failures
                try:
                    v = decode(v)
                except Exception:
                    warnings.warn(
                        f'Failed to decode "{name}" metadata '
                        "string. Check `char_encoding` parameter."
                    )

            try:
                # For convenience, if the array contains only one single
                # entry, return this entry itself.
                v = v.item()
            except ValueError:
                v = np.squeeze(v)
            ret[name] = v
        return ret


def roi_array_to_dict(a: np.ndarray) -> List[Dict[str, List[int]]]:
    """Convert the `ROIs` structured arrays to :py:class:`dict`

    Parameters
    ----------
    a
        Structured array containing ROI data

    Returns
    -------
    One dict per ROI. Keys are "top_left", "bottom_right", and "bin",
    values are tuples whose first element is the x axis value and the
    second element is the y axis value.
    """

    dict_list = []
    a = a[["startx", "starty", "endx", "endy", "groupx", "groupy"]]
    for sx, sy, ex, ey, gx, gy in a:
        roi_dict = {
            "top_left": [int(sx), int(sy)],
            "bottom_right": [int(ex), int(ey)],
            "bin": [int(gx), int(gy)],
        }
        dict_list.append(roi_dict)
    return dict_list
