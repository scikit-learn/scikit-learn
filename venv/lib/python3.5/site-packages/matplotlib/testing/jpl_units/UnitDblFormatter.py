#===========================================================================
#
# UnitDblFormatter
#
#===========================================================================


"""UnitDblFormatter module containing class UnitDblFormatter."""

#===========================================================================
# Place all imports after here.
#
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import matplotlib.ticker as ticker
#
# Place all imports before here.
#===========================================================================

__all__ = [ 'UnitDblFormatter' ]

#===========================================================================
class UnitDblFormatter( ticker.ScalarFormatter ):
   """The formatter for UnitDbl data types.  This allows for formatting
      with the unit string.
   """
   def __init__( self, *args, **kwargs ):
      'The arguments are identical to matplotlib.ticker.ScalarFormatter.'
      ticker.ScalarFormatter.__init__( self, *args, **kwargs )

   def __call__( self, x, pos = None ):
      'Return the format for tick val x at position pos'
      if len(self.locs) == 0:
         return ''
      else:
         return '{:.12}'.format(x)

   def format_data_short( self, value ):
      "Return the value formatted in 'short' format."
      return '{:.12}'.format(value)

   def format_data( self, value ):
      "Return the value formatted into a string."
      return '{:.12}'.format(value)
