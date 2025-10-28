"""Reading and saving of images and videos."""

import warnings

from .manage_plugins import *
from .manage_plugins import _hide_plugin_deprecation_warnings
from .sift import *
from .collection import *

from ._io import *
from ._image_stack import *


with _hide_plugin_deprecation_warnings():
    reset_plugins()


__all__ = [
    "concatenate_images",
    "imread",
    "imread_collection",
    "imread_collection_wrapper",
    "imsave",
    "load_sift",
    "load_surf",
    "pop",
    "push",
    "ImageCollection",
    "MultiImage",
]


def __getattr__(name):
    if name == "available_plugins":
        warnings.warn(
            "`available_plugins` is deprecated since version 0.25 and will "
            "be removed in version 0.27. Instead, use `imageio` or other "
            "I/O packages directly.",
            category=FutureWarning,
            stacklevel=2,
        )
        return globals()["_available_plugins"]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
