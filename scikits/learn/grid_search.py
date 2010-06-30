import numpy as np

from joblib import Parallel, delayed

try:
    from itertools import product
except:
    def product(*args, **kwds):
        pools = map(tuple, args) * kwds.get('repeat', 1)
        result = [[]]
        for pool in pools:
            result = [x+[y] for x in result for y in pool]
        for prod in result:
            yield tuple(prod)

def grid(**kwargs):
    """ Generators on the combination of the various parameter lists given.

        Parameters
        -----------
        kwargs: keyword arguments, lists
            Each keyword argument must be a list of values that should
            be explored.

        Returns
        --------
        params: dictionary
            Dictionnary with the input parameters taking the various
            values succesively.

        Examples
        ---------
        >>> list(grid(a=[1, 2], b=[True, False]))
        [{'a': 1, 'b': True}, {'a': 1, 'b': False}, {'a': 2, 'b': True}, {'a': 2, 'b': False}]
    """
    keys = kwargs.keys()
    for v in product(*kwargs.values()):
        params = dict(zip(keys,v))
        yield params

def fit_grid_point(X, y, clf_factory, clf_params, cross_val_factory,
                                        loss_func, **fit_params):
    """Run fit on one set of parameters
    Returns the score and the instance of the classifier
    """
    n_samples, n_features = X.shape
    clf = clf_factory(**clf_params)
    cv = cross_val_factory(n_samples)
    y_pred = np.zeros_like(y)
    for train, test in cv:
        clf.fit(X[train], y[train], **fit_params)
        y_pred[test] = np.asarray(clf.predict(X[test])).astype(np.int)

    score = loss_func(y, y_pred)
    return clf, score

class GridSearch(object):
    """
    Object to run a grid search on the parameters of a classifier.

    Important memmbers are fit, predict.

    GridSearch implements a "fit" method and a "predict" method like
    any classifier except that the parameters of the classifier
    used to predict is optimized by cross-validation

    Parameters
    ---------
    clf_factory : object type that implements the "fit" and "predict" methods
        A object of that type is instanciated for each grid point

    params : dict
        a dictionary of parameters that are used the generate the grid

    cross_val_factory : a generator to run crossvalidation

    loss_func : function that takes 2 arguments and compares them in
        order to evaluate the performance of prediciton (small is good)

    n_jobs : int
        number of jobs to run in parallel (default 1)

    Optional Parameters
    -------------------

    Members
    -------

    Examples
    --------
    >>> import numpy as np
    >>> from scikits.learn.cross_val import LeaveOneOut
    >>> from scikits.learn.svm import SVC
    >>> X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    >>> y = np.array([1, 1, 2, 2])
    >>> parameters = {'kernel':('linear', 'rbf'), 'C':[1, 10]}
    >>> def loss_func(y1, y2):
    ...     return np.mean(y1 != y2)
    >>> clf = GridSearch(SVC, parameters, LeaveOneOut, loss_func, n_jobs=1)
    >>> print clf.fit(X, y).predict([[-0.8, -1]])
    [ 1.]
    """
    def __init__(self, clf_factory, params, cross_val_factory, loss_func,
                        fit_params={}, n_jobs=1):
        self.clf_factory = clf_factory
        self.params = params
        self.cross_val_factory = cross_val_factory
        self.loss_func = loss_func
        self.n_jobs = n_jobs
        self.fit_params = fit_params

    def fit(self, X, y, **kw):
        """Run fit with all sets of parameters
        Returns the best classifier
        """
        best_score = np.inf

        self.learner = None
        self.predict = None

        g = grid(**self.params)
        out = Parallel(n_jobs=self.n_jobs)(
            delayed(fit_grid_point)(X, y, self.clf_factory, clf_params,
                    self.cross_val_factory,
                    self.loss_func, **self.fit_params) for clf_params in g)

        for clf, score in out:
            if score < best_score:
                best_score = score
                self.learner = clf
                self.predict = clf.predict

        return self.learner

if __name__ == '__main__':

    import numpy as np
    from scikits.learn.cross_val import LeaveOneOut
    from scikits.learn.svm import SVC
    X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    y = np.array([1, 1, 2, 2])
    parameters = {'kernel':('linear', 'rbf'), 'C':[1, 10]}
    def loss_func(y1, y2):
        return np.mean(y1 != y2)
    clf = GridSearch(SVC, parameters, LeaveOneOut, loss_func, n_jobs=2)
    print clf.fit(X, y).predict([[-0.8, -1]])
