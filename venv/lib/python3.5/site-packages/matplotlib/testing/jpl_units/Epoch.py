#===========================================================================
#
# Epoch
#
#===========================================================================


"""Epoch module."""

#===========================================================================
# Place all imports after here.
#
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import math
import datetime as DT
from matplotlib.dates import date2num
#
# Place all imports before here.
#===========================================================================

#===========================================================================
class Epoch(object):
   # Frame conversion offsets in seconds
   # t(TO) = t(FROM) + allowed[ FROM ][ TO ]
   allowed = {
      "ET" : {
         "UTC" : +64.1839,
         },
      "UTC" : {
         "ET" : -64.1839,
         },
      }

   #-----------------------------------------------------------------------
   def __init__( self, frame, sec=None, jd=None, daynum=None, dt=None ):
      """Create a new Epoch object.

      Build an epoch 1 of 2 ways:

      Using seconds past a Julian date:
      #   Epoch( 'ET', sec=1e8, jd=2451545 )

      or using a matplotlib day number
      #   Epoch( 'ET', daynum=730119.5 )


      = ERROR CONDITIONS
      - If the input units are not in the allowed list, an error is thrown.

      = INPUT VARIABLES
      - frame    The frame of the epoch.  Must be 'ET' or 'UTC'
      - sec      The number of seconds past the input JD.
      - jd       The Julian date of the epoch.
      - daynum   The matplotlib day number of the epoch.
      - dt       A python datetime instance.
      """
      if ( ( sec is None and jd is not None ) or
           ( sec is not None and jd is None ) or
           ( daynum is not None and ( sec is not None or jd is not None ) ) or
           ( daynum is None and dt is None and ( sec is None or jd is None ) ) or
           ( daynum is not None and dt is not None ) or
           ( dt is not None and ( sec is not None or jd is not None ) ) or
           ( (dt is not None) and not isinstance(dt, DT.datetime) ) ):
         msg = "Invalid inputs.  Must enter sec and jd together, " \
               "daynum by itself, or dt (must be a python datetime).\n" \
               "Sec = %s\nJD  = %s\ndnum= %s\ndt  = %s" \
               % ( str( sec ), str( jd ), str( daynum ), str( dt ) )
         raise ValueError( msg )

      if frame not in self.allowed:
         msg = "Input frame '%s' is not one of the supported frames of %s" \
               % ( frame, str( list(six.iterkeys(self.allowed) ) ) )
         raise ValueError(msg)

      self._frame = frame

      if dt is not None:
         daynum = date2num( dt )

      if daynum is not None:
         # 1-JAN-0001 in JD = 1721425.5
         jd = float( daynum ) + 1721425.5
         self._jd = math.floor( jd )
         self._seconds = ( jd - self._jd ) * 86400.0

      else:
         self._seconds = float( sec )
         self._jd = float( jd )

         # Resolve seconds down to [ 0, 86400 )
         deltaDays = int( math.floor( self._seconds / 86400.0 ) )
         self._jd += deltaDays
         self._seconds -= deltaDays * 86400.0

   #-----------------------------------------------------------------------
   def convert( self, frame ):
      if self._frame == frame:
         return self

      offset = self.allowed[ self._frame ][ frame ]

      return Epoch( frame, self._seconds + offset, self._jd )

   #-----------------------------------------------------------------------
   def frame( self ):
      return self._frame

   #-----------------------------------------------------------------------
   def julianDate( self, frame ):
      t = self
      if frame != self._frame:
         t = self.convert( frame )

      return t._jd + t._seconds / 86400.0

   #-----------------------------------------------------------------------
   def secondsPast( self, frame, jd ):
      t = self
      if frame != self._frame:
         t = self.convert( frame )

      delta = t._jd - jd
      return t._seconds + delta * 86400

   #-----------------------------------------------------------------------
   def __cmp__( self, rhs ):
      """Compare two Epoch's.

      = INPUT VARIABLES
      - rhs    The Epoch to compare against.

      = RETURN VALUE
      - Returns -1 if self < rhs, 0 if self == rhs, +1 if self > rhs.
      """
      t = self
      if self._frame != rhs._frame:
         t = self.convert( rhs._frame )

      if t._jd != rhs._jd:
         return cmp( t._jd, rhs._jd )

      return cmp( t._seconds, rhs._seconds )

   #-----------------------------------------------------------------------
   def __add__( self, rhs ):
      """Add a duration to an Epoch.

      = INPUT VARIABLES
      - rhs    The Epoch to subtract.

      = RETURN VALUE
      - Returns the difference of ourselves and the input Epoch.
      """
      t = self
      if self._frame != rhs.frame():
         t = self.convert( rhs._frame )

      sec = t._seconds + rhs.seconds()

      return Epoch( t._frame, sec, t._jd )

   #-----------------------------------------------------------------------
   def __sub__( self, rhs ):
      """Subtract two Epoch's or a Duration from an Epoch.

      Valid:
      Duration = Epoch - Epoch
      Epoch = Epoch - Duration

      = INPUT VARIABLES
      - rhs    The Epoch to subtract.

      = RETURN VALUE
      - Returns either the duration between to Epoch's or the a new
        Epoch that is the result of subtracting a duration from an epoch.
      """
      # Delay-load due to circular dependencies.
      import matplotlib.testing.jpl_units as U

      # Handle Epoch - Duration
      if isinstance( rhs, U.Duration ):
         return self + -rhs

      t = self
      if self._frame != rhs._frame:
         t = self.convert( rhs._frame )

      days = t._jd - rhs._jd
      sec = t._seconds - rhs._seconds

      return U.Duration( rhs._frame, days*86400 + sec )

   #-----------------------------------------------------------------------
   def __str__( self ):
      """Print the Epoch."""
      return "%22.15e %s" % ( self.julianDate( self._frame ), self._frame )

   #-----------------------------------------------------------------------
   def __repr__( self ):
      """Print the Epoch."""
      return str( self )

   #-----------------------------------------------------------------------
   def range( start, stop, step ):
      """Generate a range of Epoch objects.

      Similar to the Python range() method.  Returns the range [
      start, stop ) at the requested step.  Each element will be a
      Epoch object.

      = INPUT VARIABLES
      - start    The starting value of the range.
      - stop     The stop value of the range.
      - step     Step to use.

      = RETURN VALUE
      - Returns a list contianing the requested Epoch values.
      """
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

#===========================================================================
