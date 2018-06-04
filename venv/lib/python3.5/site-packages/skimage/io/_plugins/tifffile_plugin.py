from ...external.tifffile import TiffFile, imsave


def imread(fname, dtype=None, **kwargs):
    """Load a tiff image from file.

    Parameters
    ----------
    fname : str or file
       File name or file-like-object.
    dtype : numpy dtype object or string specifier
       Specifies data type of array elements (Not currently used).
    kwargs : keyword pairs, optional
        Additional keyword arguments to pass through (see ``tifffile``'s
        ``imread`` function).

    Notes
    -----
    Provided by Christophe Golhke's tifffile.py [1]_, and supports many
    advanced image types including multi-page and floating point.

    References
    ----------
    .. [1] http://www.lfd.uci.edu/~gohlke/code/tifffile.py

    """
    if 'img_num' in kwargs:
        kwargs['key'] = kwargs.pop('img_num')
    with open(fname, 'rb') as f:
        tif = TiffFile(f)
        return tif.asarray(**kwargs)
