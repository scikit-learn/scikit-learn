# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from time import time

import numpy as np
import torch as xp
from tqdm import tqdm

from sklearn import config_context
from sklearn.linear_model import LogisticRegression

n_samples, n_features, n_classes = 100000, 1000, 50
device = "cuda"
n_iter = 10

X_np = np.random.rand(n_samples, n_features)
y_np = np.random.randint(0, 10, n_samples)
numpy_fit_times = []
numpy_predict_times = []
for _ in tqdm(range(n_iter), desc="Numpy"):
    lr = LogisticRegression(C=0.8, solver="lbfgs", max_iter=200)
    start = time()
    lr.fit(X_np, y_np)
    numpy_fit_times.append(round(time() - start, 3))
    start = time()
    pred = lr.predict_proba(X_np)
    numpy_predict_times.append(round(time() - start, 3))

avg_numpy_fit = round(sum(numpy_fit_times) / n_iter, 3)
avg_numpy_predict = round(sum(numpy_predict_times) / n_iter, 3)

torch_fit_times = []
torch_predict_times = []
X_xp = xp.rand((n_samples, n_features), device=device)
y_xp = xp.randint(0, n_classes, (n_samples,), device=device)
for _ in tqdm(range(n_iter), desc=f"Torch {device}"):
    with config_context(array_api_dispatch=True):
        lr = LogisticRegression(C=0.8, solver="lbfgs", max_iter=200)
        start = time()
        lr.fit(X_xp, y_xp)
        torch_fit_times.append(round(time() - start, 3))
        start = time()
        pred = lr.predict_proba(X_xp)
        first = float(pred[0, 0])
        torch_predict_times.append(round(time() - start, 3))

avg_torch_fit = round(sum(torch_fit_times) / n_iter, 3)
avg_torch_predict = round(sum(torch_predict_times) / n_iter, 3)

print(f"Average fit time numpy: {avg_numpy_fit}")
print(f"Average fit time torch {device}: {avg_torch_fit}")
print(f"Torch {device} fit speedup: {round(avg_numpy_fit / avg_torch_fit, 2)}X")


print(f"Average predict time numpy: {avg_numpy_predict}")
print(f"Average predict time torch {device}: {avg_torch_predict}")
print(
    f"Torch {device} predict speedup: {round(avg_numpy_predict / avg_torch_predict, 2)}"
    "X"
)
