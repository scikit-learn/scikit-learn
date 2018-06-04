from __future__ import absolute_import, division, print_function

from functools import partial

from . import converters, exceptions, filters, validators
from ._config import get_run_validators, set_run_validators
from ._funcs import asdict, assoc, astuple, evolve, has
from ._make import (
    NOTHING, Attribute, Factory, attrib, attrs, fields, fields_dict,
    make_class, validate
)


__version__ = "18.1.0"

__title__ = "attrs"
__description__ = "Classes Without Boilerplate"
__uri__ = "http://www.attrs.org/"
__doc__ = __description__ + " <" + __uri__ + ">"

__author__ = "Hynek Schlawack"
__email__ = "hs@ox.cx"

__license__ = "MIT"
__copyright__ = "Copyright (c) 2015 Hynek Schlawack"


s = attributes = attrs
ib = attr = attrib
dataclass = partial(attrs, auto_attribs=True)  # happy Easter ;)

__all__ = [
    "Attribute",
    "Factory",
    "NOTHING",
    "asdict",
    "assoc",
    "astuple",
    "attr",
    "attrib",
    "attributes",
    "attrs",
    "converters",
    "evolve",
    "exceptions",
    "fields",
    "fields_dict",
    "filters",
    "get_run_validators",
    "has",
    "ib",
    "make_class",
    "s",
    "set_run_validators",
    "validate",
    "validators",
]
