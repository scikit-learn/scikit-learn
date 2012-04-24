import numpy as np
cimport numpy as np

cdef class LossFunction:
    """Base class for loss functions"""

    cpdef np.ndarray loss(self, np.ndarray p, np.ndarray y):
        """Evaluate the loss function.

        Parameters
        ----------
        p : np.ndarray, shape = [n_samples]
            The prediction.
        y : np.ndarray, shape = [n_samples]
            The true value (aka target).

        Returns
        -------

        np.ndarray, shape = [n_samples]
            The loss evaluated at `p` and `y`.
        """
        raise NotImplementedError()

    cpdef np.ndarray dloss(self, np.ndarray p, np.ndarray y):
        """Evaluate the derivative of the loss function with respect to
        the prediction `p`.

        Parameters
        ----------
        p : np.ndarray, shape = [n_samples]
            The prediction.
        y : np.ndarray, shape = [n_samples]
            The true value (aka target).

        Returns
        -------

        np.ndarray, shape = [n_samples]
            The derivative of the loss function at `p` and `y`.
        """
        raise NotImplementedError()


cdef class SquaredLoss(LossFunction):
    """Squared loss function."""

    cpdef np.ndarray loss(self, np.ndarray p, np.ndarray y):
        return 0.5 * (p - y) * (p - y)

    cpdef np.ndarray dloss(self, np.ndarray p, np.ndarray y):
        return p - y

    def __reduce__(self):
        return SquaredLoss, ()
