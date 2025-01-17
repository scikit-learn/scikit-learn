import functools
import warnings
from itertools import product

import numpy as np

from .dtype import img_as_float


def _rename_image_params(func):
    wm_images = (
        "Since version 0.24, the two input images are named `image0` and "
        "`image1` (instead of `image1` and `image2`, respectively). Please use "
        "`image0, image1` to avoid this warning for now, and avoid an error "
        "from version 0.26 onwards."
    )

    wm_method = (
        "Starting in version 0.24, all arguments following `image0, image1` "
        "(including `method`) will be keyword-only. Please pass `method=` "
        "in the function call to avoid this warning for now, and avoid an error "
        "from version 0.26 onwards."
    )

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # Turn all args into kwargs
        for i, (value, param) in enumerate(
            zip(args, ["image0", "image1", "method", "n_tiles"])
        ):
            if i >= 2:
                warnings.warn(wm_method, category=FutureWarning, stacklevel=2)
            if param in kwargs:
                raise ValueError(
                    f"{param} passed both as positional and keyword argument."
                )
            else:
                kwargs[param] = value
        args = tuple()

        # Account for `image2` if given
        if "image2" in kwargs.keys():
            warnings.warn(wm_images, category=FutureWarning, stacklevel=2)

            # Safely move `image2` to `image1` if that's empty
            if "image1" in kwargs.keys():
                # Safely move `image1` to `image0`
                if "image0" in kwargs.keys():
                    raise ValueError(
                        "Three input images given; please use only `image0` "
                        "and `image1`."
                    )
                kwargs["image0"] = kwargs.pop("image1")
            kwargs["image1"] = kwargs.pop("image2")

        return func(*args, **kwargs)

    return wrapper


@_rename_image_params
def compare_images(image0, image1, *, method='diff', n_tiles=(8, 8)):
    """
    Return an image showing the differences between two images.

    .. versionadded:: 0.16

    Parameters
    ----------
    image0, image1 : ndarray, shape (M, N)
        Images to process, must be of the same shape.

        .. versionchanged:: 0.24
            `image1` and `image2` were renamed into `image0` and `image1`
            respectively.
    method : string, optional
        Method used for the comparison.
        Valid values are {'diff', 'blend', 'checkerboard'}.
        Details are provided in the note section.

        .. versionchanged:: 0.24
            This parameter and following ones are keyword-only.
    n_tiles : tuple, optional
        Used only for the `checkerboard` method. Specifies the number
        of tiles (row, column) to divide the image.

    Returns
    -------
    comparison : ndarray, shape (M, N)
        Image showing the differences.

    Notes
    -----
    ``'diff'`` computes the absolute difference between the two images.
    ``'blend'`` computes the mean value.
    ``'checkerboard'`` makes tiles of dimension `n_tiles` that display
    alternatively the first and the second image. Note that images must be
    2-dimensional to be compared with the checkerboard method.
    """

    if image1.shape != image0.shape:
        raise ValueError('Images must have the same shape.')

    img1 = img_as_float(image0)
    img2 = img_as_float(image1)

    if method == 'diff':
        comparison = np.abs(img2 - img1)
    elif method == 'blend':
        comparison = 0.5 * (img2 + img1)
    elif method == 'checkerboard':
        if img1.ndim != 2:
            raise ValueError(
                'Images must be 2-dimensional to be compared with the '
                'checkerboard method.'
            )
        shapex, shapey = img1.shape
        mask = np.full((shapex, shapey), False)
        stepx = int(shapex / n_tiles[0])
        stepy = int(shapey / n_tiles[1])
        for i, j in product(range(n_tiles[0]), range(n_tiles[1])):
            if (i + j) % 2 == 0:
                mask[i * stepx : (i + 1) * stepx, j * stepy : (j + 1) * stepy] = True
        comparison = np.zeros_like(img1)
        comparison[mask] = img1[mask]
        comparison[~mask] = img2[~mask]
    else:
        raise ValueError(
            'Wrong value for `method`. '
            'Must be either "diff", "blend" or "checkerboard".'
        )
    return comparison
