"""
The classes here provide support for using custom classes with
Matplotlib, e.g., those that do not expose the array interface but know
how to convert themselves to arrays.  It also supports classes with
units and units conversion.  Use cases include converters for custom
objects, e.g., a list of datetime objects, as well as for objects that
are unit aware.  We don't assume any particular units implementation;
rather a units implementation must provide the register with the Registry
converter dictionary and a `ConversionInterface`.  For example,
here is a complete implementation which supports plotting with native
datetime objects::

    import matplotlib.units as units
    import matplotlib.dates as dates
    import matplotlib.ticker as ticker
    import datetime

    class DateConverter(units.ConversionInterface):

        @staticmethod
        def convert(value, unit, axis):
            'Convert a datetime value to a scalar or array'
            return dates.date2num(value)

        @staticmethod
        def axisinfo(unit, axis):
            'Return major and minor tick locators and formatters'
            if unit!='date': return None
            majloc = dates.AutoDateLocator()
            majfmt = dates.AutoDateFormatter(majloc)
            return AxisInfo(majloc=majloc,
                            majfmt=majfmt,
                            label='date')

        @staticmethod
        def default_units(x, axis):
            'Return the default unit for x or None'
            return 'date'

    # Finally we register our object type with the Matplotlib units registry.
    units.registry[datetime.date] = DateConverter()

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


import six
from matplotlib.cbook import iterable, is_numlike, safe_first_element
import numpy as np


class AxisInfo(object):
    """
    Information to support default axis labeling, tick labeling, and
    default limits. An instance of this class must be returned by
    :meth:`ConversionInterface.axisinfo`.
    """
    def __init__(self, majloc=None, minloc=None,
                 majfmt=None, minfmt=None, label=None,
                 default_limits=None):
        """
        Parameters
        ----------
        majloc, minloc : Locator, optional
            Tick locators for the major and minor ticks.
        majfmt, minfmt : Formatter, optional
            Tick formatters for the major and minor ticks.
        label : str, optional
            The default axis label.
        default_limits : optional
            The default min and max limits of the axis if no data has
            been plotted.

        Notes
        -----
        If any of the above are ``None``, the axis will simply use the
        default value.
        """
        self.majloc = majloc
        self.minloc = minloc
        self.majfmt = majfmt
        self.minfmt = minfmt
        self.label = label
        self.default_limits = default_limits


class ConversionInterface(object):
    """
    The minimal interface for a converter to take custom data types (or
    sequences) and convert them to values Matplotlib can use.
    """
    @staticmethod
    def axisinfo(unit, axis):
        """
        Return an `~units.AxisInfo` instance for the axis with the
        specified units.
        """
        return None

    @staticmethod
    def default_units(x, axis):
        """
        Return the default unit for *x* or ``None`` for the given axis.
        """
        return None

    @staticmethod
    def convert(obj, unit, axis):
        """
        Convert *obj* using *unit* for the specified *axis*.
        If *obj* is a sequence, return the converted sequence.
        The output must be a sequence of scalars that can be used by the numpy
        array layer.
        """
        return obj

    @staticmethod
    def is_numlike(x):
        """
        The Matplotlib datalim, autoscaling, locators etc work with
        scalars which are the units converted to floats given the
        current unit.  The converter may be passed these floats, or
        arrays of them, even when units are set.
        """
        if iterable(x):
            for thisx in x:
                return is_numlike(thisx)
        else:
            return is_numlike(x)


class Registry(dict):
    """
    A register that maps types to conversion interfaces.
    """
    def __init__(self):
        dict.__init__(self)
        self._cached = {}

    def get_converter(self, x):
        """
        Get the converter for data that has the same type as *x*. If no
        converters are registered for *x*, returns ``None``.
        """

        if not len(self):
            return None  # nothing registered
        # DISABLED idx = id(x)
        # DISABLED cached = self._cached.get(idx)
        # DISABLED if cached is not None: return cached

        converter = None
        classx = getattr(x, '__class__', None)

        if classx is not None:
            converter = self.get(classx)

        if converter is None and hasattr(x, "values"):
            # this unpacks pandas series or dataframes...
            x = x.values

        # If x is an array, look inside the array for data with units
        if isinstance(x, np.ndarray) and x.size:
            xravel = x.ravel()
            try:
                # pass the first value of x that is not masked back to
                # get_converter
                if not np.all(xravel.mask):
                    # some elements are not masked
                    converter = self.get_converter(
                        xravel[np.argmin(xravel.mask)])
                    return converter
            except AttributeError:
                # not a masked_array
                # Make sure we don't recurse forever -- it's possible for
                # ndarray subclasses to continue to return subclasses and
                # not ever return a non-subclass for a single element.
                next_item = xravel[0]
                if (not isinstance(next_item, np.ndarray) or
                        next_item.shape != x.shape):
                    converter = self.get_converter(next_item)
                return converter

        # If we haven't found a converter yet, try to get the first element
        if converter is None:
            try:
                thisx = safe_first_element(x)
            except (TypeError, StopIteration):
                pass
            else:
                if classx and classx != getattr(thisx, '__class__', None):
                    converter = self.get_converter(thisx)
                    return converter

        # DISABLED self._cached[idx] = converter
        return converter


registry = Registry()
