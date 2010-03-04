import numpy as np
from math import exp,log
import scipy

def fast_logdet(A):
    """
    Compute log(det(A)) for A symmetric
    Equivalent to : np.log(nl.det(A))
    but more robust
    It returns -Inf if det(A) is non positive or is not defined.
    """
    ld = np.sum(np.log(np.diag(A)))
    a = exp(ld/A.shape[0])
    d = scipy.linalg.det(A/a)
    ld += log(d)
    if not np.isfinite(ld):
        return -np.inf
    return ld
