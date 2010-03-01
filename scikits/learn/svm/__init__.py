"""
Support Vector Machine module.


Examples
--------
>>> from scikits.learn import svm
>>> import numpy as np
>>> data = np.array([[-1, -1], [-2, -2], [1, 1], [2, 2]])
>>> labels = np.array([1, 1, 2, 2])
>>> print svm.predict(data, labels, [[1, 0]])
[ 2.]

See also
--------
:func:`scikits.learn.svm.predict`
:class:`scikits.learn.svm.SVM`
"""
from svm import SVM, predict
