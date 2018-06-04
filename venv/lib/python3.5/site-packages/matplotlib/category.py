# -*- coding: utf-8 -*-
"""
Module that allows plotting of string "category" data.  i.e.
``plot(['d', 'f', 'a'],[1, 2, 3])`` will plot three points with x-axis
values of 'd', 'f', 'a'.

See :doc:`/gallery/lines_bars_and_markers/categorical_variables` for an
example.

The module uses Matplotlib's `matplotlib.units` mechanism to convert from
strings to integers, provides a tick locator and formatter, and the
class:`.UnitData` that creates and stores the string-to-integer mapping.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from collections import OrderedDict
import itertools

import six


import numpy as np

import matplotlib.units as units
import matplotlib.ticker as ticker

# np 1.6/1.7 support
from distutils.version import LooseVersion

VALID_TYPES = tuple(set(six.string_types +
                        (bytes, six.text_type, np.str_, np.bytes_)))


class StrCategoryConverter(units.ConversionInterface):
    @staticmethod
    def convert(value, unit, axis):
        """Converts strings in value to floats using
        mapping information store in the  unit object

        Parameters
        ----------
        value : string or iterable
            value or list of values to be converted
        unit : :class:`.UnitData`
           object string unit information for value
        axis : :class:`~matplotlib.Axis.axis`
            axis on which the converted value is plotted

        Returns
        -------
        mapped_ value : float or ndarray[float]

        .. note:: axis is not used in this function
        """
        # dtype = object preserves numerical pass throughs
        values = np.atleast_1d(np.array(value, dtype=object))

        # pass through sequence of non binary numbers
        if all((units.ConversionInterface.is_numlike(v) and
                not isinstance(v, VALID_TYPES)) for v in values):
            return np.asarray(values, dtype=float)

        # force an update so it also does type checking
        unit.update(values)

        str2idx = np.vectorize(unit._mapping.__getitem__,
                               otypes=[float])

        mapped_value = str2idx(values)
        return mapped_value

    @staticmethod
    def axisinfo(unit, axis):
        """Sets the default axis ticks and labels

        Parameters
        ---------
        unit : :class:`.UnitData`
            object string unit information for value
        axis : :class:`~matplotlib.Axis.axis`
            axis for which information is being set

        Returns
        -------
        :class:~matplotlib.units.AxisInfo~
            Information to support default tick labeling

        .. note: axis is not used
        """
        # locator and formatter take mapping dict because
        # args need to be pass by reference for updates
        majloc = StrCategoryLocator(unit._mapping)
        majfmt = StrCategoryFormatter(unit._mapping)
        return units.AxisInfo(majloc=majloc, majfmt=majfmt)

    @staticmethod
    def default_units(data, axis):
        """ Sets and updates the :class:`~matplotlib.Axis.axis~ units

        Parameters
        ----------
        data : string or iterable of strings
        axis : :class:`~matplotlib.Axis.axis`
            axis on which the data is plotted

        Returns
        -------
        class:~.UnitData~
            object storing string to integer mapping
        """
        # the conversion call stack is supposed to be
        # default_units->axis_info->convert
        if axis.units is None:
            axis.set_units(UnitData(data))
        else:
            axis.units.update(data)
        return axis.units


class StrCategoryLocator(ticker.Locator):
    """tick at every integer mapping of the string data"""
    def __init__(self, units_mapping):
        """
        Parameters
        -----------
        units_mapping : Dict[str, int]
             string:integer mapping
        """
        self._units = units_mapping

    def __call__(self):
        return list(self._units.values())

    def tick_values(self, vmin, vmax):
        return self()


class StrCategoryFormatter(ticker.Formatter):
    """String representation of the data at every tick"""
    def __init__(self, units_mapping):
        """
        Parameters
        ----------
        units_mapping : Dict[Str, int]
            string:integer mapping
        """
        self._units = units_mapping

    def __call__(self, x, pos=None):
        if pos is None:
            return ""
        r_mapping = {v: StrCategoryFormatter._text(k)
                     for k, v in self._units.items()}
        return r_mapping.get(int(np.round(x)), '')

    @staticmethod
    def _text(value):
        """Converts text values into `utf-8` or `ascii` strings
        """
        if LooseVersion(np.__version__) < LooseVersion('1.7.0'):
            if (isinstance(value, (six.text_type, np.unicode))):
                value = value.encode('utf-8', 'ignore').decode('utf-8')
        if isinstance(value, (np.bytes_, six.binary_type)):
            value = value.decode(encoding='utf-8')
        elif not isinstance(value, (np.str_, six.string_types)):
            value = str(value)
        return value


class UnitData(object):
    def __init__(self, data=None):
        """Create mapping between unique categorical values
        and integer identifiers
        ----------
        data: iterable
              sequence of string values
        """
        self._mapping = OrderedDict()
        self._counter = itertools.count(start=0)
        if data is not None:
            self.update(data)

    def update(self, data):
        """Maps new values to integer identifiers.

        Paramters
        ---------
        data: iterable
              sequence of string values

        Raises
        ------
        TypeError
              If the value in data is not a string, unicode, bytes type
        """
        data = np.atleast_1d(np.array(data, dtype=object))

        for val in OrderedDict.fromkeys(data):
            if not isinstance(val, VALID_TYPES):
                raise TypeError("{val!r} is not a string".format(val=val))
            if val not in self._mapping:
                self._mapping[val] = next(self._counter)


# Connects the convertor to matplotlib
units.registry[str] = StrCategoryConverter()
units.registry[np.str_] = StrCategoryConverter()
units.registry[six.text_type] = StrCategoryConverter()
units.registry[bytes] = StrCategoryConverter()
units.registry[np.bytes_] = StrCategoryConverter()
