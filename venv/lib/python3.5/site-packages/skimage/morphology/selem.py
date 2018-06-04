import numpy as np
from scipy import ndimage as ndi
from .. import draw


def square(width, dtype=np.uint8):
    """Generates a flat, square-shaped structuring element.

    Every pixel along the perimeter has a chessboard distance
    no greater than radius (radius=floor(width/2)) pixels.

    Parameters
    ----------
    width : int
        The width and height of the square.

    Other Parameters
    ----------------
    dtype : data-type
        The data type of the structuring element.

    Returns
    -------
    selem : ndarray
        A structuring element consisting only of ones, i.e. every
        pixel belongs to the neighborhood.

    """
    return np.ones((width, width), dtype=dtype)


def rectangle(width, height, dtype=np.uint8):
    """Generates a flat, rectangular-shaped structuring element.

    Every pixel in the rectangle generated for a given width and given height
    belongs to the neighborhood.

    Parameters
    ----------
    width : int
        The width of the rectangle.
    height : int
        The height of the rectangle.

    Other Parameters
    ----------------
    dtype : data-type
        The data type of the structuring element.

    Returns
    -------
    selem : ndarray
        A structuring element consisting only of ones, i.e. every
        pixel belongs to the neighborhood.

    """
    return np.ones((width, height), dtype=dtype)


def diamond(radius, dtype=np.uint8):
    """Generates a flat, diamond-shaped structuring element.

    A pixel is part of the neighborhood (i.e. labeled 1) if
    the city block/Manhattan distance between it and the center of
    the neighborhood is no greater than radius.

    Parameters
    ----------
    radius : int
        The radius of the diamond-shaped structuring element.

    Other Parameters
    ----------------
    dtype : data-type
        The data type of the structuring element.

    Returns
    -------

    selem : ndarray
        The structuring element where elements of the neighborhood
        are 1 and 0 otherwise.
    """
    L = np.arange(0, radius * 2 + 1)
    I, J = np.meshgrid(L, L)
    return np.array(np.abs(I - radius) + np.abs(J - radius) <= radius,
                    dtype=dtype)


def disk(radius, dtype=np.uint8):
    """Generates a flat, disk-shaped structuring element.

    A pixel is within the neighborhood if the euclidean distance between
    it and the origin is no greater than radius.

    Parameters
    ----------
    radius : int
        The radius of the disk-shaped structuring element.

    Other Parameters
    ----------------
    dtype : data-type
        The data type of the structuring element.

    Returns
    -------
    selem : ndarray
        The structuring element where elements of the neighborhood
        are 1 and 0 otherwise.
    """
    L = np.arange(-radius, radius + 1)
    X, Y = np.meshgrid(L, L)
    return np.array((X ** 2 + Y ** 2) <= radius ** 2, dtype=dtype)


def ellipse(width, height, dtype=np.uint8):
    """Generates a flat, ellipse-shaped structuring element.

    Every pixel along the perimeter of ellipse satisfies
    the equation ``(x/width+1)**2 + (y/height+1)**2 = 1``.

    Parameters
    ----------
    width : int
        The width of the ellipse-shaped structuring element.
    height : int
        The height of the ellipse-shaped structuring element.

    Other Parameters
    ----------------
    dtype : data-type
        The data type of the structuring element.

    Returns
    -------
    selem : ndarray
        The structuring element where elements of the neighborhood
        are 1 and 0 otherwise.

    Examples
    --------
    >>> from skimage.morphology import selem
    >>> selem.ellipse(5, 3)
    array([[0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0],
           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
           [0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0]], dtype=uint8)

    """
    selem = np.zeros((2 * height + 1, 2 * width + 1), dtype=dtype)
    rows, cols = draw.ellipse(height, width, height + 1, width + 1)
    selem[rows, cols] = 1
    return selem


def cube(width, dtype=np.uint8):
    """ Generates a cube-shaped structuring element.

    This is the 3D equivalent of a square.
    Every pixel along the perimeter has a chessboard distance
    no greater than radius (radius=floor(width/2)) pixels.

    Parameters
    ----------
    width : int
        The width, height and depth of the cube.

    Other Parameters
    ----------------
    dtype : data-type
        The data type of the structuring element.

    Returns
    -------
    selem : ndarray
        A structuring element consisting only of ones, i.e. every
        pixel belongs to the neighborhood.

    """
    return np.ones((width, width, width), dtype=dtype)


def octahedron(radius, dtype=np.uint8):
    """Generates a octahedron-shaped structuring element.

    This is the 3D equivalent of a diamond.
    A pixel is part of the neighborhood (i.e. labeled 1) if
    the city block/Manhattan distance between it and the center of
    the neighborhood is no greater than radius.

    Parameters
    ----------
    radius : int
        The radius of the octahedron-shaped structuring element.

    Other Parameters
    ----------------
    dtype : data-type
        The data type of the structuring element.

    Returns
    -------

    selem : ndarray
        The structuring element where elements of the neighborhood
        are 1 and 0 otherwise.
    """
    # note that in contrast to diamond(), this method allows non-integer radii
    n = 2 * radius + 1
    Z, Y, X = np.mgrid[-radius:radius:n * 1j,
                       -radius:radius:n * 1j,
                       -radius:radius:n * 1j]
    s = np.abs(X) + np.abs(Y) + np.abs(Z)
    return np.array(s <= radius, dtype=dtype)


def ball(radius, dtype=np.uint8):
    """Generates a ball-shaped structuring element.

    This is the 3D equivalent of a disk.
    A pixel is within the neighborhood if the euclidean distance between
    it and the origin is no greater than radius.

    Parameters
    ----------
    radius : int
        The radius of the ball-shaped structuring element.

    Other Parameters
    ----------------
    dtype : data-type
        The data type of the structuring element.

    Returns
    -------
    selem : ndarray
        The structuring element where elements of the neighborhood
        are 1 and 0 otherwise.
    """
    n = 2 * radius + 1
    Z, Y, X = np.mgrid[-radius:radius:n * 1j,
                       -radius:radius:n * 1j,
                       -radius:radius:n * 1j]
    s = X ** 2 + Y ** 2 + Z ** 2
    return np.array(s <= radius * radius, dtype=dtype)


def octagon(m, n, dtype=np.uint8):
    """Generates an octagon shaped structuring element.

    For a given size of (m) horizontal and vertical sides
    and a given (n) height or width of slanted sides octagon is generated.
    The slanted sides are 45 or 135 degrees to the horizontal axis
    and hence the widths and heights are equal.

    Parameters
    ----------
    m : int
        The size of the horizontal and vertical sides.
    n : int
        The height or width of the slanted sides.

    Other Parameters
    ----------------
    dtype : data-type
        The data type of the structuring element.

    Returns
    -------
    selem : ndarray
        The structuring element where elements of the neighborhood
        are 1 and 0 otherwise.

    """
    from . import convex_hull_image
    selem = np.zeros((m + 2 * n, m + 2 * n))
    selem[0, n] = 1
    selem[n, 0] = 1
    selem[0, m + n - 1] = 1
    selem[m + n - 1, 0] = 1
    selem[-1, n] = 1
    selem[n, -1] = 1
    selem[-1, m + n - 1] = 1
    selem[m + n - 1, -1] = 1
    selem = convex_hull_image(selem).astype(dtype)
    return selem


def star(a, dtype=np.uint8):
    """Generates a star shaped structuring element.

    Start has 8 vertices and is an overlap of square of size `2*a + 1`
    with its 45 degree rotated version.
    The slanted sides are 45 or 135 degrees to the horizontal axis.

    Parameters
    ----------
    a : int
        Parameter deciding the size of the star structural element. The side
        of the square array returned is `2*a + 1 + 2*floor(a / 2)`.

    Other Parameters
    ----------------
    dtype : data-type
        The data type of the structuring element.

    Returns
    -------
    selem : ndarray
        The structuring element where elements of the neighborhood
        are 1 and 0 otherwise.

    """
    from . import convex_hull_image

    if a == 1:
        bfilter = np.zeros((3, 3), dtype)
        bfilter[:] = 1
        return bfilter

    m = 2 * a + 1
    n = a // 2
    selem_square = np.zeros((m + 2 * n, m + 2 * n))
    selem_square[n: m + n, n: m + n] = 1

    c = (m + 2 * n - 1) // 2
    selem_rotated = np.zeros((m + 2 * n, m + 2 * n))
    selem_rotated[0, c] = selem_rotated[-1, c] = 1
    selem_rotated[c, 0] = selem_rotated[c, -1] = 1
    selem_rotated = convex_hull_image(selem_rotated).astype(int)

    selem = selem_square + selem_rotated
    selem[selem > 0] = 1

    return selem.astype(dtype)


def _default_selem(ndim):
    """Generates a cross-shaped structuring element (connectivity=1).

    This is the default structuring element (selem) if no selem was specified.

    Parameters
    ----------
    ndim : int
        Number of dimensions of the image.

    Returns
    -------
    selem : ndarray
        The structuring element where elements of the neighborhood
        are 1 and 0 otherwise.

    """
    return ndi.morphology.generate_binary_structure(ndim, 1)
