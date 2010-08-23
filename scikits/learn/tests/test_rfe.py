"""
Testing Recursive feature elimination

"""

import numpy as np
from scikits.learn.svm import SVC
from scikits.learn.cross_val import StratifiedKFold
from scikits.learn import datasets
from scikits.learn.rfe import RFECV
from scikits.learn.metrics import zero_one

################################################################################
# Loading a dataset
iris = datasets.load_iris()
X = iris.data
y = iris.target

# Some noisy data not correlated
random = np.random.RandomState(seed=0)
E = random.normal(size=(len(X), 5))

# Add the noisy data to the informative features
X = np.c_[X, E]

def test_rfe():
    """Check that rfe recoverse the correct features on IRIS dataset"""

    svc = SVC(kernel='linear')
    rfecv = RFECV(estimator=svc, n_features=4, percentage=0.1, loss_func=zero_one)
    rfecv.fit(X, y, cv=StratifiedKFold(y, 2))
    X_r = rfecv.transform(X)

    assert rfecv.support_.sum() == iris.data.shape[1]

