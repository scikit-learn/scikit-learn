from __future__ import absolute_import, division, print_function

import numpy as np
from toolz import merge, partial

from ..base import tokenize
from .. import threaded


def _partial_fit(model, x, y, kwargs=None):
    kwargs = kwargs or dict()
    model.partial_fit(x, y, **kwargs)
    return model


def fit(model, x, y, get=threaded.get, **kwargs):
    """ Fit scikit learn model against dask arrays

    Model must support the ``partial_fit`` interface for online or batch
    learning.

    This method will be called on dask arrays in sequential order.  Ideally
    your rows are independent and identically distributed.

    Parameters
    ----------
    model: sklearn model
        Any model supporting partial_fit interface
    x: dask Array
        Two dimensional array, likely tall and skinny
    y: dask Array
        One dimensional array with same chunks as x's rows
    kwargs:
        options to pass to partial_fit

    Examples
    --------
    >>> import dask.array as da
    >>> X = da.random.random((10, 3), chunks=(5, 3))
    >>> y = da.random.randint(0, 2, 10, chunks=(5,))

    >>> from sklearn.linear_model import SGDClassifier
    >>> sgd = SGDClassifier()

    >>> sgd = da.learn.fit(sgd, X, y, classes=[1, 0])
    >>> sgd  # doctest: +SKIP
    SGDClassifier(alpha=0.0001, class_weight=None, epsilon=0.1, eta0=0.0,
           fit_intercept=True, l1_ratio=0.15, learning_rate='optimal',
           loss='hinge', n_iter=5, n_jobs=1, penalty='l2', power_t=0.5,
           random_state=None, shuffle=False, verbose=0, warm_start=False)

    This passes all of X and y through the classifier sequentially.  We can use
    the classifier as normal on in-memory data

    >>> import numpy as np
    >>> sgd.predict(np.random.random((4, 3)))  # doctest: +SKIP
    array([1, 0, 0, 1])

    Or predict on a larger dataset

    >>> z = da.random.random((400, 3), chunks=(100, 3))
    >>> da.learn.predict(sgd, z)  # doctest: +SKIP
    dask.array<x_11, shape=(400,), chunks=((100, 100, 100, 100),), dtype=int64>
    """
    assert x.ndim == 2
    if y is not None:
        assert y.ndim == 1
        assert x.chunks[0] == y.chunks[0]
    assert hasattr(model, 'partial_fit')
    if len(x.chunks[1]) > 1:
        x = x.reblock(chunks=(x.chunks[0], sum(x.chunks[1])))

    nblocks = len(x.chunks[0])

    name = 'fit-' + tokenize(model, x, y, kwargs)
    dsk = {(name, -1): model}
    dsk.update(dict(((name, i), (_partial_fit, (name, i - 1),
                                               (x.name, i, 0),
                                               (getattr(y, 'name', ''), i), kwargs))
                    for i in range(nblocks)))

    return get(merge(x.dask, getattr(y, 'dask', {}), dsk), (name, nblocks - 1))


def _predict(model, x):
    return model.predict(x)[:, None]


def predict(model, x):
    """ Predict with a scikit learn model

    Parameters
    ----------
    model : scikit learn classifier
    x : dask Array

    See docstring for ``da.learn.fit``
    """
    assert x.ndim == 2
    if len(x.chunks[1]) > 1:
        x = x.reblock(chunks=(x.chunks[0], sum(x.chunks[1])))
    func = partial(_predict, model)
    xx = np.zeros((1, x.shape[1]), dtype=x.dtype)
    dt = model.predict(xx).dtype
    return x.map_blocks(func, chunks=(x.chunks[0], (1,)), dtype=dt).squeeze()
