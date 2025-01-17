# -*- coding: utf-8 -*-
# imageio is distributed under the terms of the (new) BSD License.

# flake8: noqa

"""
Here you can find documentation on how to write your own plugin to allow
ImageIO to access a new backend. Plugins are quite object oriented, and
the relevant classes and their interaction are documented here:

.. currentmodule:: imageio

.. autosummary::
    :toctree: ../_autosummary
    :template: better_class.rst

    imageio.core.Format
    imageio.core.Request

.. note::
    You can always check existing plugins if you want to see examples.

What methods to implement
-------------------------

To implement a new plugin, create a new class that inherits from
:class:`imageio.core.Format`. and implement the following functions:

.. autosummary::
    :toctree: ../_autosummary

    imageio.core.Format.__init__
    imageio.core.Format._can_read
    imageio.core.Format._can_write

Further, each format contains up to two nested classes; one for reading and
one for writing. To support reading and/or writing, the respective classes
need to be defined.

For reading, create a nested class that inherits from
``imageio.core.Format.Reader`` and that implements the following functions:

    * Implement ``_open(**kwargs)`` to initialize the reader. Deal with the
        user-provided keyword arguments here.
    * Implement ``_close()`` to clean up.
    * Implement ``_get_length()`` to provide a suitable length based on what
        the user expects. Can be ``inf`` for streaming data.
    * Implement ``_get_data(index)`` to return an array and a meta-data dict.
    * Implement ``_get_meta_data(index)`` to return a meta-data dict. If index
        is None, it should return the 'global' meta-data.

For writing, create a nested class that inherits from
``imageio.core.Format.Writer`` and implement the following functions:

    * Implement ``_open(**kwargs)`` to initialize the writer. Deal with the
        user-provided keyword arguments here.
    * Implement ``_close()`` to clean up.
    * Implement ``_append_data(im, meta)`` to add data (and meta-data).
    * Implement ``_set_meta_data(meta)`` to set the global meta-data.

"""

import importlib
import os
import warnings


# v2 imports remove in v3
from .. import formats

# v2 allows formatting plugins by environment variable
# this is done here.
env_plugin_order = os.getenv("IMAGEIO_FORMAT_ORDER", None)
if env_plugin_order is not None:  # pragma: no cover
    warnings.warn(
        "Setting plugin priority through an environment variable is"
        " deprecated and will be removed in ImageIO v3. There is no"
        " replacement planned for this feature. If you have an"
        " active use-case for it, please reach out to us on GitHub.",
        DeprecationWarning,
    )

    formats.sort(*os.getenv("IMAGEIO_FORMAT_ORDER", "").split(","))


# this class replaces plugin module. For details
# see https://stackoverflow.com/questions/2447353/getattr-on-a-module
def __getattr__(name):
    """Lazy-Import Plugins

    This function dynamically loads plugins into the imageio.plugin
    namespace upon first access. For example, the following snippet will
    delay importing freeimage until the second line:

    >>> import imageio
    >>> imageio.plugins.freeimage.download()

    """

    try:
        return importlib.import_module(f"imageio.plugins.{name}")
    except ImportError:
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'") from None
