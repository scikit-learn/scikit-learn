from time import time

from sklearn.datasets import make_regression, make_classification
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.gbm import GBMRegressor
from sklearn.gbm import GBMClassifier

import pstats
import cProfile
import pygbm

classif = True
n_classes = 2
n_samples = int(1e6)
max_iter = 5

if classif:
    X, y = make_classification(n_samples=n_samples, random_state=0, n_classes=n_classes, n_clusters_per_class=1)
    GBM = GBMClassifier
    GBDT = GradientBoostingClassifier
    PYGBM_GBM = pygbm.GradientBoostingClassifier
else:
    X, y = make_regression(n_samples=n_samples, random_state=0)
    GBM = GBMRegressor
    GBDT = GradientBoostingRegressor
    PYGBM_GBM = pygbm.GradientBoostingRegressor


pygbm_est = PYGBM_GBM(
    max_iter=max_iter,
    scoring=None,  # no early stopping
    validation_split=None,
    random_state=0,
    verbose=False)
print("compiling pygbm code")
pygbm_est.fit(X[:1000], y[:1000])
print("done")

gbm = GBM(
    max_iter=max_iter,
    scoring=None,  # no early stopping
    validation_split=None,
    n_iter_no_change=None,
    random_state=0,
    verbose=True)
tic = time()
gbm.fit(X, y)
fit_duration = time() - tic
tic = time()
print(f'score: {gbm.score(X, y)}')
score_duration = time() - tic
print(f'sklearn gbm fit_duration: {fit_duration:.3f}s')
print(f'sklearn gbm score_duration {score_duration:.3f}s')


pygbm_est.set_params(verbose=True)
tic = time()
pygbm_est.fit(X, y)
fit_duration = time() - tic
tic = time()
print(f'score: {pygbm_est.score(X, y)}')
score_duration = time() - tic
print(f'pygbm fit_duration: {fit_duration:.3f}s')
print(f'pygbm score_duration {score_duration:.3f}s')

# cProfile.runctx("gbm.fit(X, y)", globals(), locals(), "Profile.prof")
# s = pstats.Stats("Profile.prof")
# s.strip_dirs().sort_stats("time").print_stats(.2)

