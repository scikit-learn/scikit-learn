# Copyright (c) 2006-2012 Filip Wasilewski <http://en.ig.ma/>
# Copyright (c) 2012-2016 The PyWavelets Developers
#                         <https://github.com/PyWavelets/pywt>
# See COPYING for license details.

"""
The thresholding helper module implements the most popular signal thresholding
functions.
"""

from __future__ import division, print_function, absolute_import

__all__ = ['threshold']

import numpy as np


def soft(data, value, substitute=0):
    data = np.asarray(data)

    magnitude = np.absolute(data)
    sign = np.sign(data)
    thresholded = (magnitude - value).clip(0) * sign

    cond = np.less(magnitude, value)
    return np.where(cond, substitute, thresholded)


def hard(data, value, substitute=0):
    data = np.asarray(data)
    cond = np.less(np.absolute(data), value)
    return np.where(cond, substitute, data)


def greater(data, value, substitute=0):
    data = np.asarray(data)
    return np.where(np.less(data, value), substitute, data)


def less(data, value, substitute=0):
    data = np.asarray(data)
    return np.where(np.greater(data, value), substitute, data)


thresholding_options = {'soft': soft,
                        'hard': hard,
                        'greater': greater,
                        'less': less}


def threshold(data, value, mode='soft', substitute=0):
    """
    Thresholds the input data depending on the mode argument.

    In ``soft`` thresholding, the data values where their absolute value is
    less than the value param are replaced with substitute. From the data
    values with absolute value greater or equal to the thresholding value,
    a quantity of ``(signum * value)`` is subtracted.

    In ``hard`` thresholding, the data values where their absolute value is
    less than the value param are replaced with substitute. Data values with
    absolute value greater or equal to the thresholding value stay untouched.

    In ``greater`` thresholding, the data is replaced with substitute where
    data is below the thresholding value. Greater data values pass untouched.

    In ``less`` thresholding, the data is replaced with substitute where data
    is above the thresholding value. Less data values pass untouched.

    Parameters
    ----------
    data : array_like
        Numeric data.
    value : scalar
        Thresholding value.
    mode : {'soft', 'hard', 'greater', 'less'}
        Decides the type of thresholding to be applied on input data. Default
        is 'soft'.
    substitute : float, optional
        Substitute value (default: 0).

    Returns
    -------
    output : array
        Thresholded array.

    Examples
    --------
    >>> import numpy as np
    >>> import pywt
    >>> data = np.linspace(1, 4, 7)
    >>> data
    array([ 1. ,  1.5,  2. ,  2.5,  3. ,  3.5,  4. ])
    >>> pywt.threshold(data, 2, 'soft')
    array([ 0. ,  0. ,  0. ,  0.5,  1. ,  1.5,  2. ])
    >>> pywt.threshold(data, 2, 'hard')
    array([ 0. ,  0. ,  2. ,  2.5,  3. ,  3.5,  4. ])
    >>> pywt.threshold(data, 2, 'greater')
    array([ 0. ,  0. ,  2. ,  2.5,  3. ,  3.5,  4. ])
    >>> pywt.threshold(data, 2, 'less')
    array([ 1. ,  1.5,  2. ,  0. ,  0. ,  0. ,  0. ])

    """

    try:
        return thresholding_options[mode](data, value, substitute)
    except KeyError:
        # Make sure error is always identical by sorting keys
        keys = ("'{0}'".format(key) for key in
                sorted(thresholding_options.keys()))
        raise ValueError("The mode parameter only takes values from: {0}."
                         .format(', '.join(keys)))
