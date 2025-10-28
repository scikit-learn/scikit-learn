from tifffile import imread as tifffile_imread
from tifffile import imwrite as tifffile_imwrite

__all__ = ['imread', 'imsave']


def imsave(fname, arr, **kwargs):
    """Load a tiff image to file.

    Parameters
    ----------
    fname : str or file
        File name or file-like object.
    arr : ndarray
        The array to write.
    kwargs : keyword pairs, optional
        Additional keyword arguments to pass through (see ``tifffile``'s
        ``imwrite`` function).

    Notes
    -----
    Provided by the tifffile library [1]_, and supports many
    advanced image types including multi-page and floating-point.

    This implementation will set ``photometric='RGB'`` when writing if the first
    or last axis of `arr` has length 3 or 4. To override this, explicitly
    pass the ``photometric`` kwarg.

    This implementation will set ``planarconfig='SEPARATE'`` when writing if the
    first axis of arr has length 3 or 4. To override this, explicitly
    specify the ``planarconfig`` kwarg.

    References
    ----------
    .. [1] https://pypi.org/project/tifffile/

    """
    if arr.shape[0] in [3, 4]:
        if 'planarconfig' not in kwargs:
            kwargs['planarconfig'] = 'SEPARATE'
        rgb = True
    else:
        rgb = arr.shape[-1] in [3, 4]
    if rgb and 'photometric' not in kwargs:
        kwargs['photometric'] = 'RGB'

    return tifffile_imwrite(fname, arr, **kwargs)


def imread(fname, **kwargs):
    """Load a tiff image from file.

    Parameters
    ----------
    fname : str or file
        File name or file-like-object.
    kwargs : keyword pairs, optional
        Additional keyword arguments to pass through (see ``tifffile``'s
        ``imread`` function).

    Notes
    -----
    Provided by the tifffile library [1]_, and supports many
    advanced image types including multi-page and floating point.

    References
    ----------
    .. [1] https://pypi.org/project/tifffile/

    """
    if 'img_num' in kwargs:
        kwargs['key'] = kwargs.pop('img_num')

    return tifffile_imread(fname, **kwargs)
