#
# The Python Imaging Library
# $Id$
#
# path interface
#
# History:
# 1996-11-04 fl   Created
# 2002-04-14 fl   Added documentation stub class
#
# Copyright (c) Secret Labs AB 1997.
# Copyright (c) Fredrik Lundh 1996.
#
# See the README file for information on usage and redistribution.
#

from . import Image


# the Python class below is overridden by the C implementation.


class Path(object):

    def __init__(self, xy):
        pass

    def compact(self, distance=2):
        """
        Compacts the path, by removing points that are close to each other.
        This method modifies the path in place.
        """
        pass

    def getbbox(self):
        """Gets the bounding box."""
        pass

    def map(self, function):
        """Maps the path through a function."""
        pass

    def tolist(self, flat=0):
        """
        Converts the path to Python list.
        #
        @param flat By default, this function returns a list of 2-tuples
            [(x, y), ...].  If this argument is true, it returns a flat list
            [x, y, ...] instead.
        @return A list of coordinates.
        """
        pass

    def transform(self, matrix):
        """Transforms the path."""
        pass


# override with C implementation
Path = Image.core.path
