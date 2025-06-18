# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from time import time

import numpy as np
import torch as xp
from tqdm import tqdm

from sklearn._config import config_context
from sklearn.preprocessing._polynomial import PolynomialFeatures

X_np = np.random.rand(100000, 100)
X_xp_cuda = xp.asarray(X_np, device="cuda")

# Numpy benchmarks
fit_times = []
transform_times = []
for _ in tqdm(range(10), desc="Numpy Flow"):
    start = time()
    pf_np = PolynomialFeatures(degree=2)
    pf_np.fit(X_np)
    fit_times.append(time() - start)

    start = time()
    pf_np.transform(X_np)
    transform_times.append(time() - start)

avg_fit_time = sum(fit_times) / 10
avg_transform_time = sum(transform_times) / 10
print(f"Avg fit time for numpy: {avg_fit_time}")
print(f"Avg transform time for numpy: {avg_transform_time}")


# Torch cuda benchmarks
fit_times = []
transform_times = []
for _ in tqdm(range(10), desc="Torch cuda Flow"):
    with config_context(array_api_dispatch=True):
        start = time()
        pf_xp = PolynomialFeatures(degree=2)
        pf_xp.fit(X_xp_cuda)
        fit_times.append(time() - start)

        start = time()
        pf_xp.transform(X_xp_cuda)
        transform_times.append(time() - start)

avg_fit_time = sum(fit_times) / 10
avg_transform_time = sum(transform_times) / 10
print(f"Avg fit time for torch cuda: {avg_fit_time}")
print(f"Avg transform time for torch cuda: {avg_transform_time}")
