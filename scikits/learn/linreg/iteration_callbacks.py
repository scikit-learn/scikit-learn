"""Iterations callbacks

    Iterations callbacks can be used for early stopping in an iterative
    algorithm or for example for vizualisation purposes during the
    optimization.

"""

class BaseIterationCallback(object):
    """Base callback to be called at the end of each iteration of CD

    To be subclassed if more monitoring is required.
    """

    def __init__(self, linear_model):
        pass

    def __call__(self, n_iter, **kwargs):
        return True


class IterationCallbackMaxit(BaseIterationCallback):
    """Callback to be called at the end of each iteration of CD

    - record and check the duality gap for early stop of the optim
      (before maxiter)

    """

    def __init__(self, maxit):
        self.maxit = maxit

    def __call__(self, n_iter, **kwargs):
        return n_iter < self.maxit


class IterationCallbackFunc(BaseIterationCallback):
    """Callback to be called at the end of each iteration of CD

    - Store the values computed during iteration process

    """

    def __init__(self, func, tol=None, record=True):
        self.func = func
        self.values = []
        self.record = record
        self.tol = tol

    def __call__(self, n_iter, **kwargs):
        if n_iter == 0:
            self.values = []
            return True

        val = self.func(**kwargs)
        if self.record:
            self.values.append(val)
        else:
            self.values = val

        if self.tol is not None:
            return abs(val) > self.tol

        return True

class IterationCallbackMaxTime(BaseIterationCallback):
    """Callback to be called at the end of each iteration of CD

    - Stop convergence after a certain number of seconds of computation

    """

    def __init__(self, tmax=None):
        self.tmax = tmax
        self.t0 = None

    def __call__(self, n_iter, **kwargs):
        if n_iter == 0:
            self.t0 = time.time()
            return True

        if (time.time() - self.t0) > self.tmax:
            return False
        else:
            return True

