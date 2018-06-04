#===========================================================================
#
# UnitDbl
#
#===========================================================================


"""UnitDbl module."""

#===========================================================================
# Place all imports after here.
#
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
#
# Place all imports before here.
#===========================================================================


#===========================================================================
class UnitDbl(object):
   """Class UnitDbl in development.
   """
   #-----------------------------------------------------------------------
   # Unit conversion table.  Small subset of the full one but enough
   # to test the required functions.  First field is a scale factor to
   # convert the input units to the units of the second field.  Only
   # units in this table are allowed.
   allowed = {
               "m" : ( 0.001, "km" ),
               "km" : ( 1, "km" ),
               "mile" : ( 1.609344, "km" ),

               "rad" : ( 1, "rad" ),
               "deg" : ( 1.745329251994330e-02, "rad" ),

               "sec" : ( 1, "sec" ),
               "min" : ( 60.0, "sec" ),
               "hour" : ( 3600, "sec" ),
             }

   _types = {
              "km" : "distance",
              "rad" : "angle",
              "sec" : "time",
            }

   #-----------------------------------------------------------------------
   def __init__( self, value, units ):
      """Create a new UnitDbl object.

      Units are internally converted to km, rad, and sec.  The only
      valid inputs for units are [ m, km, mile, rad, deg, sec, min, hour ].

      The field UnitDbl.value will contain the converted value.  Use
      the convert() method to get a specific type of units back.

      = ERROR CONDITIONS
      - If the input units are not in the allowed list, an error is thrown.

      = INPUT VARIABLES
      - value    The numeric value of the UnitDbl.
      - units    The string name of the units the value is in.
      """
      self.checkUnits( units )

      data = self.allowed[ units ]
      self._value = float( value * data[0] )
      self._units = data[1]

   #-----------------------------------------------------------------------
   def convert( self, units ):
      """Convert the UnitDbl to a specific set of units.

      = ERROR CONDITIONS
      - If the input units are not in the allowed list, an error is thrown.

      = INPUT VARIABLES
      - units    The string name of the units to convert to.

      = RETURN VALUE
      - Returns the value of the UnitDbl in the requested units as a floating
        point number.
      """
      if self._units == units:
         return self._value

      self.checkUnits( units )

      data = self.allowed[ units ]
      if self._units != data[1]:
         msg = "Error trying to convert to different units.\n" \
               "   Invalid conversion requested.\n" \
               "   UnitDbl: %s\n" \
               "   Units:   %s\n" % ( str( self ), units )
         raise ValueError( msg )

      return self._value / data[0]

   #-----------------------------------------------------------------------
   def __abs__( self ):
      """Return the absolute value of this UnitDbl."""
      return UnitDbl( abs( self._value ), self._units )

   #-----------------------------------------------------------------------
   def __neg__( self ):
      """Return the negative value of this UnitDbl."""
      return UnitDbl( -self._value, self._units )

   #-----------------------------------------------------------------------
   def __nonzero__( self ):
      """Test a UnitDbl for a non-zero value.

      = RETURN VALUE
      - Returns true if the value is non-zero.
      """
      if six.PY3:
          return self._value.__bool__()
      else:
          return self._value.__nonzero__()

   if six.PY3:
      __bool__ = __nonzero__

   #-----------------------------------------------------------------------
   def __cmp__( self, rhs ):
      """Compare two UnitDbl's.

      = ERROR CONDITIONS
      - If the input rhs units are not the same as our units,
        an error is thrown.

      = INPUT VARIABLES
      - rhs    The UnitDbl to compare against.

      = RETURN VALUE
      - Returns -1 if self < rhs, 0 if self == rhs, +1 if self > rhs.
      """
      self.checkSameUnits( rhs, "compare" )
      return cmp( self._value, rhs._value )

   #-----------------------------------------------------------------------
   def __add__( self, rhs ):
      """Add two UnitDbl's.

      = ERROR CONDITIONS
      - If the input rhs units are not the same as our units,
        an error is thrown.

      = INPUT VARIABLES
      - rhs    The UnitDbl to add.

      = RETURN VALUE
      - Returns the sum of ourselves and the input UnitDbl.
      """
      self.checkSameUnits( rhs, "add" )
      return UnitDbl( self._value + rhs._value, self._units )

   #-----------------------------------------------------------------------
   def __sub__( self, rhs ):
      """Subtract two UnitDbl's.

      = ERROR CONDITIONS
      - If the input rhs units are not the same as our units,
        an error is thrown.

      = INPUT VARIABLES
      - rhs    The UnitDbl to subtract.

      = RETURN VALUE
      - Returns the difference of ourselves and the input UnitDbl.
      """
      self.checkSameUnits( rhs, "subtract" )
      return UnitDbl( self._value - rhs._value, self._units )

   #-----------------------------------------------------------------------
   def __mul__( self, rhs ):
      """Scale a UnitDbl by a value.

      = INPUT VARIABLES
      - rhs    The scalar to multiply by.

      = RETURN VALUE
      - Returns the scaled UnitDbl.
      """
      return UnitDbl( self._value * rhs, self._units )

   #-----------------------------------------------------------------------
   def __rmul__( self, lhs ):
      """Scale a UnitDbl by a value.

      = INPUT VARIABLES
      - lhs    The scalar to multiply by.

      = RETURN VALUE
      - Returns the scaled UnitDbl.
      """
      return UnitDbl( self._value * lhs, self._units )

   #-----------------------------------------------------------------------
   def __div__( self, rhs ):
      """Divide a UnitDbl by a value.

      = INPUT VARIABLES
      - rhs    The scalar to divide by.

      = RETURN VALUE
      - Returns the scaled UnitDbl.
      """
      return UnitDbl( self._value / rhs, self._units )

   #-----------------------------------------------------------------------
   def __str__( self ):
      """Print the UnitDbl."""
      return "%g *%s" % ( self._value, self._units )

   #-----------------------------------------------------------------------
   def __repr__( self ):
      """Print the UnitDbl."""
      return "UnitDbl( %g, '%s' )" % ( self._value, self._units )

   #-----------------------------------------------------------------------
   def type( self ):
      """Return the type of UnitDbl data."""
      return self._types[ self._units ]

   #-----------------------------------------------------------------------
   def range( start, stop, step=None ):
      """Generate a range of UnitDbl objects.

      Similar to the Python range() method.  Returns the range [
      start, stop ) at the requested step.  Each element will be a
      UnitDbl object.

      = INPUT VARIABLES
      - start    The starting value of the range.
      - stop     The stop value of the range.
      - step     Optional step to use.  If set to None, then a UnitDbl of
                 value 1 w/ the units of the start is used.

      = RETURN VALUE
      - Returns a list contianing the requested UnitDbl values.
      """
      if step is None:
         step = UnitDbl( 1, start._units )

      elems = []

      i = 0
      while True:
         d = start + i * step
         if d >= stop:
            break

         elems.append( d )
         i += 1

      return elems

   range = staticmethod( range )

   #-----------------------------------------------------------------------
   def checkUnits( self, units ):
      """Check to see if some units are valid.

      = ERROR CONDITIONS
      - If the input units are not in the allowed list, an error is thrown.

      = INPUT VARIABLES
      - units    The string name of the units to check.
      """
      if units not in self.allowed:
         msg = "Input units '%s' are not one of the supported types of %s" \
               % ( units, str( list(six.iterkeys(self.allowed)) ) )
         raise ValueError( msg )

   #-----------------------------------------------------------------------
   def checkSameUnits( self, rhs, func ):
      """Check to see if units are the same.

      = ERROR CONDITIONS
      - If the units of the rhs UnitDbl are not the same as our units,
        an error is thrown.

      = INPUT VARIABLES
      - rhs    The UnitDbl to check for the same units
      - func   The name of the function doing the check.
      """
      if self._units != rhs._units:
         msg = "Cannot %s units of different types.\n" \
               "LHS: %s\n" \
               "RHS: %s" % ( func, self._units, rhs._units )
         raise ValueError( msg )

#===========================================================================
