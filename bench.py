# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from time import time

import numpy as np
import torch
from tqdm import tqdm

from sklearn.base import config_context
from sklearn.calibration import _TemperatureScaling

xp = torch
device_ = "cuda"
dtype_np = np.float64
dtype_xp = xp.float64
n_samples = 1000000
n_classes = 300


execution_times = []
for _ in tqdm(range(10), desc="Numpy"):
    y = np.random.randint(0, n_classes, n_samples).astype(dtype_np)
    pred = np.random.rand(n_samples, n_classes).astype(dtype_np)
    start = time()
    cal = _TemperatureScaling()
    cal.fit(pred, y)
    execution_times.append(time() - start)

avg_time = sum(execution_times) / 10
print(f"Avg execution_time for numpy: {avg_time}")

execution_times = []
for _ in tqdm(range(10), desc=f"Torch {device_}"):
    y = xp.randint(0, n_classes, (n_samples,), dtype=dtype_xp, device=device_)
    pred = xp.rand((n_samples, n_classes), dtype=dtype_xp, device=device_)
    start = time()
    with config_context(array_api_dispatch=True):
        cal = _TemperatureScaling()
        cal.fit(pred, y)
    execution_times.append(time() - start)

avg_time = sum(execution_times) / 10
print(f"Avg execution_time for torch {device_}: {avg_time}")
