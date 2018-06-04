#===========================================================================
#
# StrConverter
#
#===========================================================================


"""StrConverter module containing class StrConverter."""

#===========================================================================
# Place all imports after here.
#
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import xrange

import matplotlib.units as units
from matplotlib.cbook import iterable

# Place all imports before here.
#===========================================================================

__all__ = [ 'StrConverter' ]

#===========================================================================
class StrConverter( units.ConversionInterface ):
   """: A matplotlib converter class.  Provides matplotlib conversion
        functionality for string data values.

   Valid units for string are:
   - 'indexed' : Values are indexed as they are specified for plotting.
   - 'sorted'  : Values are sorted alphanumerically.
   - 'inverted' : Values are inverted so that the first value is on top.
   - 'sorted-inverted' :  A combination of 'sorted' and 'inverted'
   """

   #------------------------------------------------------------------------
   @staticmethod
   def axisinfo( unit, axis ):
      """: Returns information on how to handle an axis that has string data.

      = INPUT VARIABLES
      - axis    The axis using this converter.
      - unit    The units to use for a axis with string data.

      = RETURN VALUE
      - Returns a matplotlib AxisInfo data structure that contains
        minor/major formatters, major/minor locators, and default
        label information.
      """

      return None

   #------------------------------------------------------------------------
   @staticmethod
   def convert( value, unit, axis ):
      """: Convert value using unit to a float.  If value is a sequence, return
      the converted sequence.

      = INPUT VARIABLES
      - axis    The axis using this converter.
      - value   The value or list of values that need to be converted.
      - unit    The units to use for a axis with Epoch data.

      = RETURN VALUE
      - Returns the value parameter converted to floats.
      """

      if ( units.ConversionInterface.is_numlike( value ) ):
         return value

      if ( value == [] ):
         return []

      # we delay loading to make matplotlib happy
      ax = axis.axes
      if axis is ax.get_xaxis():
         isXAxis = True
      else:
         isXAxis = False

      axis.get_major_ticks()
      ticks = axis.get_ticklocs()
      labels = axis.get_ticklabels()

      labels = [ l.get_text() for l in labels if l.get_text() ]

      if ( not labels ):
         ticks = []
         labels = []


      if ( not iterable( value ) ):
         value = [ value ]

      newValues = []
      for v in value:
         if ( (v not in labels) and (v not in newValues) ):
            newValues.append( v )

      for v in newValues:
         if ( labels ):
            labels.append( v )
         else:
            labels = [ v ]

      #DISABLED: This is disabled because matplotlib bar plots do not
      #DISABLED: recalculate the unit conversion of the data values
      #DISABLED: this is due to design and is not really a bug.
      #DISABLED: If this gets changed, then we can activate the following
      #DISABLED: block of code.  Note that this works for line plots.
      #DISABLED if ( unit ):
      #DISABLED    if ( unit.find( "sorted" ) > -1 ):
      #DISABLED       labels.sort()
      #DISABLED    if ( unit.find( "inverted" ) > -1 ):
      #DISABLED       labels = labels[ ::-1 ]

      # add padding (so they do not appear on the axes themselves)
      labels = [ '' ] + labels + [ '' ]
      ticks = list(xrange( len(labels) ))
      ticks[0] = 0.5
      ticks[-1] = ticks[-1] - 0.5

      axis.set_ticks( ticks )
      axis.set_ticklabels( labels )
      # we have to do the following lines to make ax.autoscale_view work
      loc = axis.get_major_locator()
      loc.set_bounds( ticks[0], ticks[-1] )

      if ( isXAxis ):
         ax.set_xlim( ticks[0], ticks[-1] )
      else:
         ax.set_ylim( ticks[0], ticks[-1] )

      result = []
      for v in value:
         # If v is not in labels then something went wrong with adding new
         # labels to the list of old labels.
         errmsg  = "This is due to a logic error in the StrConverter class.  "
         errmsg += "Please report this error and its message in bugzilla."
         assert ( v in labels ), errmsg
         result.append( ticks[ labels.index(v) ] )

      ax.viewLim.ignore(-1)
      return result

   #------------------------------------------------------------------------
   @staticmethod
   def default_units( value, axis ):
      """: Return the default unit for value, or None.

      = INPUT VARIABLES
      - axis    The axis using this converter.
      - value   The value or list of values that need units.

      = RETURN VALUE
      - Returns the default units to use for value.
      Return the default unit for value, or None.
      """

      # The default behavior for string indexing.
      return "indexed"
