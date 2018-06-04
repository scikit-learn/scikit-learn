"""
skimage.novice
==============
A special Python image submodule for beginners.

Description
-----------
``skimage.novice`` provides a simple image manipulation interface for
beginners.  It allows for easy loading, manipulating, and saving of image
files.

This module is primarily intended for teaching and differs significantly from
the normal, array-oriented image functions used by scikit-image.

.. note::

    This module uses the Cartesian coordinate system, where the origin is at
    the lower-left corner instead of the upper-right and the order is x, y
    instead of row, column.


Examples
--------
We can create a Picture object open opening an image file:

>>> from skimage import novice
>>> from skimage import data
>>> picture = novice.open(data.data_dir + '/chelsea.png')

We can display pictures (after running this command, close the window to access the prompt again):

>>> picture.show()  # doctest: +SKIP

Pictures know their format:

>>> picture.format
'png'

... and where they came from:

>>> picture.path.endswith('chelsea.png')
True

... and their size:

>>> picture.size
(451, 300)
>>> picture.width
451

As a reminder, we can preview the picture with our earlier command:

>>> picture.show()  # doctest: +SKIP

Changing `size` resizes the picture.

>>> picture.size = (45, 30)

We can preview the changes we made to the picture with the ``compare`` command:

>>> picture.compare()  # doctest: +SKIP

You can iterate over pixels, which have RGB values between 0 and 255,
and know their location in the picture.

>>> for pixel in picture:
...     if (pixel.red > 128) and (pixel.x < picture.width):
...         pixel.red = pixel.red / 2

Pictures know if they've been modified from the original file

>>> picture.modified
True
>>> print(picture.path)
None

Pictures can be indexed like arrays

>>> picture[0:20, 0:20] = (0, 0, 0)

Saving the picture updates the path attribute, format, and modified state.

>>> picture.save('save-demo.jpg')
>>> picture.path.endswith('save-demo.jpg')
True
>>> picture.format
'jpeg'
>>> picture.modified
False

An image can also be restored to its original state after modification:

>>> picture[0:20, 0:20] = (0, 0, 0)
>>> picture.compare()  # doctest: +SKIP
>>> picture.reset()
>>> picture.compare()  # doctest: +SKIP

"""
import warnings
from ._novice import Picture, open, colors, color_dict


warnings.warn("The `skimage.novice` module was deprecated in version 0.14. "
              "It will be removed in 0.16.")

__all__ = ['Picture', 'open', 'colors', 'color_dict']
