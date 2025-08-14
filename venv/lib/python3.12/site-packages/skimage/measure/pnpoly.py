from ._pnpoly import _grid_points_in_poly, _points_in_poly


def grid_points_in_poly(shape, verts, binarize=True):
    """Test whether points on a specified grid are inside a polygon.

    For each ``(r, c)`` coordinate on a grid, i.e. ``(0, 0)``, ``(0, 1)`` etc.,
    test whether that point lies inside a polygon.

    You can control the output type with the `binarize` flag. Please refer to its
    documentation for further details.

    Parameters
    ----------
    shape : tuple (M, N)
        Shape of the grid.
    verts : (V, 2) array
        Specify the V vertices of the polygon, sorted either clockwise
        or anti-clockwise. The first point may (but does not need to be)
        duplicated.
    binarize: bool
        If `True`, the output of the function is a boolean mask.
        Otherwise, it is a labeled array. The labels are:
        O - outside, 1 - inside, 2 - vertex, 3 - edge.

    See Also
    --------
    points_in_poly

    Returns
    -------
    mask : (M, N) ndarray
        If `binarize` is True, the output is a boolean mask. True means the
        corresponding pixel falls inside the polygon.
        If `binarize` is False, the output is a labeled array, with pixels
        having a label between 0 and 3. The meaning of the values is:
        O - outside, 1 - inside, 2 - vertex, 3 - edge.

    """
    output = _grid_points_in_poly(shape, verts)
    if binarize:
        output = output.astype(bool)
    return output


def points_in_poly(points, verts):
    """Test whether points lie inside a polygon.

    Parameters
    ----------
    points : (K, 2) array
        Input points, ``(x, y)``.
    verts : (L, 2) array
        Vertices of the polygon, sorted either clockwise or anti-clockwise.
        The first point may (but does not need to be) duplicated.

    See Also
    --------
    grid_points_in_poly

    Returns
    -------
    mask : (K,) array of bool
        True if corresponding point is inside the polygon.

    """
    return _points_in_poly(points, verts)
