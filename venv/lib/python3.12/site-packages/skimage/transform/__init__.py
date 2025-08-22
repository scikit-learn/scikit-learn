"""Geometric and other transformations, e.g., rotations, Radon transform.

- Geometric transformation:
  These transforms change the shape or position of an image.
  They are useful for tasks such as image registration,
  alignment, and geometric correction.
  Examples: :class:`~skimage.transform.AffineTransform`,
  :class:`~skimage.transform.ProjectiveTransform`,
  :class:`~skimage.transform.EuclideanTransform`.

- Image resizing and rescaling:
  These transforms change the size or resolution of an image.
  They are useful for tasks such as down-sampling an image to
  reduce its size or up-sampling an image to increase its resolution.
  Examples: :func:`~skimage.transform.resize`,
  :func:`~skimage.transform.rescale`.

- Feature detection and extraction:
  These transforms identify and extract specific features or
  patterns in an image. They are useful for tasks such as object
  detection, image segmentation, and  feature matching.
  Examples: :func:`~skimage.transform.hough_circle`,
  :func:`~skimage.transform.pyramid_expand`,
  :func:`~skimage.transform.radon`.

- Image transformation:
  These transforms change the appearance of an image without changing its
  content. They are useful for tasks such a creating image mosaics,
  applying artistic effects, and visualizing image data.
  Examples: :func:`~skimage.transform.warp`,
  :func:`~skimage.transform.iradon`.

"""

import lazy_loader as _lazy

__getattr__, __dir__, __all__ = _lazy.attach_stub(__name__, __file__)
