from time import time

from sklearn.datasets import make_regression, make_classification
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import GBMRegressor
from sklearn.ensemble import GBMClassifier

import pstats
import cProfile

classif = False
n_samples = 500000
max_iter = 5

if classif:
    X, y = make_classification(n_samples=n_samples, random_state=0)
    GBM = GBMClassifier
    GBDT = GradientBoostingClassifier
else:
    X, y = make_regression(n_samples=n_samples, random_state=0)
    GBM = GBMRegressor
    GBDT = GradientBoostingRegressor


tic = time()
gbm = GBM(max_iter=max_iter,
                   scoring=None,  # no early stopping
                   validation_split=None,
                   n_iter_no_change=None,
                   random_state=0,
                   verbose=True)
gbm.fit(X, y)
print(f'score: {gbm.score(X, y)}')
duration = time() - tic
print(f'Took {duration:.3f}s\n')

# cProfile.runctx("gbm.fit(X, y).predict(X)", globals(), locals(), "Profile.prof")

# s = pstats.Stats("Profile.prof")
# s.strip_dirs().sort_stats("time").print_stats(.2)

tic = time()
gbdt = GBDT(n_estimators=max_iter,
            n_iter_no_change=None,  # no early stopping
            random_state=0,
            verbose=True).fit(X, y)
print(gbdt.n_estimators_)
print(f'score: {gbdt.score(X, y)}')
duration = time() - tic
print(f'Took {duration:.3f}s')
