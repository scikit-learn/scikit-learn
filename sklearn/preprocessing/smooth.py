# Authors: Luigi De Bianchi <luigi.debianchi@gmail.com>
# License: BSD 3 clause

from __future__ import division

import numpy as np

from ..base import BaseEstimator, TransformerMixin
from ..utils import check_array
from ..utils.validation import check_is_fitted

POSSIBLE_MODES = ['mirror', 'constant', 'nearest', 'wrap', 'interp']

__all__ = [
    'RectangularSmoother',
    'rectangular_smoother',
]


def rectangular_smoother(x, size=3, mode='constant', cval=0.0):
    """Smooth a features set using the "rectangular" or "unweighted
    sliding-average smooth" method.

    Parameters
    ----------
    X : {array-like}, shape [n_samples, n_features]
        The data to smooth per feature.

    size : any odd number, optional, default 3
        The window size the algorithm use to determine the new value of a point.

    mode : 'mirror', 'constant', 'nearest', 'wrap', 'interp', optional,
            default constant
        This determines the type of extension to use for the padded signal to
        which the filter is applied. When mode is ‘constant’, the padding value
        is given by cval. When the ‘interp’ mode is selected, no extension is
        used.
        Details on the mode options:
            ‘mirror’:
                Repeats the values at the edges in reverse order. The value
                closest to the edge is not included.
            ‘nearest’:
                The extension contains the nearest input value.
            ‘constant’:
                The extension contains the value given by the cval argument.
            ‘wrap’:
                The extension contains the values from the other end of the
                array.
        E.g.
            mode       |   Ext   |         Input          |   Ext
            -----------+---------+------------------------+---------
            'mirror'   | 4  3  2 | 1  2  3  4  5  6  7  8 | 7  6  5
            'nearest'  | 1  1  1 | 1  2  3  4  5  6  7  8 | 8  8  8
            'constant' | 0  0  0 | 1  2  3  4  5  6  7  8 | 0  0  0
            'wrap'     | 6  7  8 | 1  2  3  4  5  6  7  8 | 1  2  3

    cval : any number, optional, default 0.0
        Only used by mode constant. Define the constant value appendend to
        array during the evaluation.

    Returns
    -------
    X : array-like, shape [n_samples, n_features]
        Smoothed input X.

    See also
    --------
    RectangularSmoother: Performs normalization using the ``Transformer`` API
        (e.g. as part of a preprocessing :class:`sklearn.pipeline.Pipeline`).

    Notes
    -----
    For a comparison of the different scalers, transformers, and normalizers,
    see :ref:`examples/preprocessing/plot_all_scaling.py
    <sphx_glr_auto_examples_preprocessing_plot_all_scaling.py>`.

    """
    if mode not in POSSIBLE_MODES:
        raise ValueError("Given mode is not supported.")

    if size % 2 == 0:
        raise ValueError("Size must be odd.")

    whisker = int((size - 1) / 2)
    x = check_array(x, ensure_min_samples=whisker, warn_on_dtype=True,
                    estimator='The rectangular_smoother function.')

    array_filler_up, array_filler_down = populate_fillers(x, mode, whisker,
                                                          cval)

    supported_input = np.concatenate((array_filler_up, x, array_filler_down),
                                     axis=0)
    if mode == 'interp':
        result = np.zeros((x.shape[0] - size + 1, x.shape[1]))
    else:
        result = np.zeros(x.shape)

    result[0, :] = sum_samples(supported_input, 0, size)
    for row in range(1, result.shape[0]):
        result[row, :] = result[row - 1, :] \
                         - supported_input[row - 1, :] \
                         + supported_input[row + size - 1, :]
    result = np.divide(result, size)
    return result


class RectangularSmoother(BaseEstimator, TransformerMixin):
    """Smooth a features set using the "rectangular" or "unweighted
    sliding-average smooth" method.

    To smooth a data set is to create an approximating function that attempts
    to capture important patterns in the data, while leaving out noise or
    other fine-scale structures/rapid phenomena.

    Parameters
    ----------
    size : any odd number, optional, default 3
        The window size the algorithm use to determine the new value of a point.

    mode : 'mirror', 'constant', 'nearest', 'wrap', 'interp', optional,
            default constant
        This determines the type of extension to use for the padded signal to
        which the filter is applied. When mode is ‘constant’, the padding value
        is given by cval. When the ‘interp’ mode is selected, no extension is
        used.

        Details on the mode options:
            ‘mirror’:
                Repeats the values at the edges in reverse order. The value
                closest to the edge is not included.
            ‘nearest’:
                The extension contains the nearest input value.
            ‘constant’:
                The extension contains the value given by the cval argument.
            ‘wrap’:
                The extension contains the values from the other end of the
                array.
        E.g.
            mode       |   Ext   |         Input          |   Ext
            -----------+---------+------------------------+---------
            'mirror'   | 4  3  2 | 1  2  3  4  5  6  7  8 | 7  6  5
            'nearest'  | 1  1  1 | 1  2  3  4  5  6  7  8 | 8  8  8
            'constant' | 0  0  0 | 1  2  3  4  5  6  7  8 | 0  0  0
            'wrap'     | 6  7  8 | 1  2  3  4  5  6  7  8 | 1  2  3

    cval : any number, optional, default 0.0
        Only used by mode constant. Define the constant value appendend to
        array during the evaluation.

    Notes
    -----
    This estimator is stateless (besides constructor parameters), the
    fit method does nothing but is useful when used in a pipeline.

    For a comparison of the different scalers, transformers, and normalizers,
    see :ref:`examples/preprocessing/plot_all_scaling.py
    <sphx_glr_auto_examples_preprocessing_plot_all_scaling.py>`.

    """
    def __init__(self, size=3, mode='constant', cval=0.0):
        if mode not in POSSIBLE_MODES:
            raise ValueError("Given mode is not supported.")

        if size % 2 == 0:
            raise ValueError("Size must be odd.")

        self.mode = mode
        self.size = size
        self.cval = cval

    def fit(self, x, y=None):
        """Do nothing and return the estimator unchanged.

        This method is just there to implement the usual API and hence
        work in pipelines.

        Parameters
        ----------
        X : array-like
        """
        x = check_array(x,
                        ensure_min_samples=int((self.size - 1) / 2),
                        warn_on_dtype=True, estimator=self)
        self._feature_indices = x.shape[1]
        return self

    def transform(self, x):
        """Smooth each feature with the parameters defined in the constructor.
        If the mode 'interp' is used the returned array is 'size' elements
        shorter than the original one.

        Parameters
        ----------
        X : array-like, shape [n_samples, n_features]
            The data to smooth, row by row.
        y : (ignored)
            .. deprecated:: 0.19
               This parameter will be removed in 0.21.
        """
        check_is_fitted(self, 'mode')

        x = check_array(x)
        if x.shape[1] != self._feature_indices:
            raise ValueError("X has different shape than during fitting."
                             " Expected %d, got %d."
                             % (self._feature_indices, x.shape[1]))

        return rectangular_smoother(x, size=self.size, mode=self.mode,
                                    cval=self.cval)


def sum_samples(x, start, end):
    """This method evaluate the sum of the rows from 'start' index to the 'end'
    index of the passed array 'x'.

    Parameters
    ----------
    X : array-like, shape [n_samples, n_features]
        Bidimensional array

    start : int
        Starting index

    end : int
        Ending index

    Returns
    -------
    X : array-like, shape [1, n_features]
        The sum column by column of the selected indexes.
    """
    result = np.zeros(x.shape[1])
    for row in x[start:end, :]:
        result += row
    return result


def populate_fillers(x, mode, whisker, cval=0.0):
    """This method evaluate the sum of the rows from 'start' index to the 'end'
    index of the passed array 'x'.

    Parameters
    ----------
     X : {array-like}, shape [n_samples, n_features]
        The data to smooth per feature.

    mode : 'mirror', 'constant', 'nearest', 'wrap', 'interp'
        This determines the type of extension to use for the padded signal to
        which the filter is applied. When mode is ‘constant’, the padding value
        is given by cval. When the ‘interp’ mode is selected, no extension is
        used.
        Details on the mode options:
            ‘mirror’:
                Repeats the values at the edges in reverse order. The value
                closest to the edge is not included.
            ‘nearest’:
                The extension contains the nearest input value.
            ‘constant’:
                The extension contains the value given by the cval argument.
            ‘wrap’:
                The extension contains the values from the other end of the
                array.
        E.g.
            mode       |   Ext   |         Input          |   Ext
            -----------+---------+------------------------+---------
            'mirror'   | 4  3  2 | 1  2  3  4  5  6  7  8 | 7  6  5
            'nearest'  | 1  1  1 | 1  2  3  4  5  6  7  8 | 8  8  8
            'constant' | 0  0  0 | 1  2  3  4  5  6  7  8 | 0  0  0
            'wrap'     | 6  7  8 | 1  2  3  4  5  6  7  8 | 1  2  3

    whisker : int
        It's the length of the generated extensions. In general is defined as
        half the size of the current filter ( int((size - 1) / 2) ).

    cval : any number, optional, default 0.0
        Only used by mode constant. Define the constant value appendend to
        array during the evaluation.

    Returns
    -------
    filler_up : array-like, shape [whisker, n_features]
        The upper extension of the given array.

    filler_down : array-like, shape [whisker, n_features]
        The lower extension of the given array.
    """
    if x.shape[0] < whisker:
        raise ValueError("Too few sample with respect to the chosen window "
                         "size")

    if mode == 'interp':
        filler = np.zeros((0, x.shape[1]))
        return filler, filler

    filler_up = np.zeros((whisker, x.shape[1]))
    filler_down = np.zeros((whisker, x.shape[1]))

    if mode == 'mirror':
        for i in range(0, whisker):
            filler_up[i, :] = x[whisker - i, :]
            filler_down[i, :] = x[- 2 - i, :]
        return filler_up, filler_down

    if mode == 'constant':
        filler_up[:, :] = cval
        return filler_up, filler_up

    if mode == 'nearest':
        filler_up[:, :] = x[0, :]
        filler_down[:, :] = x[-1, :]
        return filler_up, filler_down

    if mode == 'wrap':
        filler_up[:, :] = x[-whisker:, :]
        filler_down[:, :] = x[:whisker, :]
        return filler_up, filler_down
