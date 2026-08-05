"""
Test performance of LogisticRegressionCV fitting with different solvers and
parallelism levels.
"""

from time import time

import numpy as np

from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import train_test_split

X, y = make_classification(
    n_samples=100_000, n_features=20, n_informative=2, n_redundant=2, random_state=42
)

train_samples = 1000  # Samples used for training the models
X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    shuffle=False,
    test_size=100_000 - train_samples,
)

print("| n_jobs | solver          | elapsed |")
print("|--------|-----------------|---------|")

for solver in [
    "lbfgs",
    "liblinear",
    "newton-cd-gram",
    "newton-cg",
    "newton-cholesky",
    "sag",
    "saga",
]:
    for n_jobs in [1, 4, 8]:
        lr = LogisticRegressionCV(
            Cs=np.logspace(-6, 6, 101),
            cv=10,
            l1_ratios=(0,),
            scoring="neg_log_loss",
            max_iter=1_000,
            use_legacy_attributes=False,
            n_jobs=n_jobs,
            solver=solver,
        )
        start = time()
        lr.fit(X_train, y_train)
        print(f"|{n_jobs:8}| {solver:16}|{round(time() - start, 3):9}|")
