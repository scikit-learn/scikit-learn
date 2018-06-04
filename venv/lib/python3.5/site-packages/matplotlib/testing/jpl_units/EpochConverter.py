#===========================================================================
#
# EpochConverter
#
#===========================================================================


"""EpochConverter module containing class EpochConverter."""

#===========================================================================
# Place all imports after here.
#
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import matplotlib.units as units
import matplotlib.dates as date_ticker
from matplotlib.cbook import iterable
#
# Place all imports before here.
#===========================================================================

__all__ = [ 'EpochConverter' ]

#===========================================================================
class EpochConverter( units.ConversionInterface ):
   """: A matplotlib converter class.  Provides matplotlib conversion
        functionality for Monte Epoch and Duration classes.
   """

   # julian date reference for "Jan 1, 0001" minus 1 day because
   # matplotlib really wants "Jan 0, 0001"
   jdRef = 1721425.5 - 1

   #------------------------------------------------------------------------
   @staticmethod
   def axisinfo( unit, axis ):
      """: Returns information on how to handle an axis that has Epoch data.

      = INPUT VARIABLES
      - unit    The units to use for a axis with Epoch data.

      = RETURN VALUE
      - Returns a matplotlib AxisInfo data structure that contains
        minor/major formatters, major/minor locators, and default
        label information.
      """

      majloc = date_ticker.AutoDateLocator()
      majfmt = date_ticker.AutoDateFormatter( majloc )

      return units.AxisInfo( majloc = majloc,
                             majfmt = majfmt,
                             label = unit )

   #------------------------------------------------------------------------
   @staticmethod
   def float2epoch( value, unit ):
      """: Convert a matplotlib floating-point date into an Epoch of the
           specified units.

      = INPUT VARIABLES
      - value    The matplotlib floating-point date.
      - unit     The unit system to use for the Epoch.

      = RETURN VALUE
      - Returns the value converted to an Epoch in the specified time system.
      """
      # Delay-load due to circular dependencies.
      import matplotlib.testing.jpl_units as U

      secPastRef = value * 86400.0 * U.UnitDbl( 1.0, 'sec' )
      return U.Epoch( unit, secPastRef, EpochConverter.jdRef )

   #------------------------------------------------------------------------
   @staticmethod
   def epoch2float( value, unit ):
      """: Convert an Epoch value to a float suitible for plotting as a
           python datetime object.

      = INPUT VARIABLES
      - value   An Epoch or list of Epochs that need to be converted.
      - unit    The units to use for an axis with Epoch data.

      = RETURN VALUE
      - Returns the value parameter converted to floats.
      """
      return value.julianDate( unit ) - EpochConverter.jdRef

   #------------------------------------------------------------------------
   @staticmethod
   def duration2float( value ):
      """: Convert a Duration value to a float suitible for plotting as a
           python datetime object.

      = INPUT VARIABLES
      - value   A Duration or list of Durations that need to be converted.

      = RETURN VALUE
      - Returns the value parameter converted to floats.
      """
      return value.seconds() / 86400.0

   #------------------------------------------------------------------------
   @staticmethod
   def convert( value, unit, axis ):
      """: Convert value using unit to a float.  If value is a sequence, return
      the converted sequence.

      = INPUT VARIABLES
      - value   The value or list of values that need to be converted.
      - unit    The units to use for an axis with Epoch data.

      = RETURN VALUE
      - Returns the value parameter converted to floats.
      """
      # Delay-load due to circular dependencies.
      import matplotlib.testing.jpl_units as U

      isNotEpoch = True
      isDuration = False

      if ( iterable(value) and not isinstance(value, six.string_types) ):
         if ( len(value) == 0 ):
            return []
         else:
            return [ EpochConverter.convert( x, unit, axis ) for x in value ]

      if ( isinstance(value, U.Epoch) ):
         isNotEpoch = False
      elif ( isinstance(value, U.Duration) ):
         isDuration = True

      if ( isNotEpoch and not isDuration and
           units.ConversionInterface.is_numlike( value ) ):
         return value

      if ( unit == None ):
         unit = EpochConverter.default_units( value, axis )

      if ( isDuration ):
         return EpochConverter.duration2float( value )
      else:
         return EpochConverter.epoch2float( value, unit )

   #------------------------------------------------------------------------
   @staticmethod
   def default_units( value, axis ):
      """: Return the default unit for value, or None.

      = INPUT VARIABLES
      - value   The value or list of values that need units.

      = RETURN VALUE
      - Returns the default units to use for value.
      """
      frame = None
      if ( iterable(value) and not isinstance(value, six.string_types) ):
         return EpochConverter.default_units( value[0], axis )
      else:
         frame = value.frame()

      return frame
