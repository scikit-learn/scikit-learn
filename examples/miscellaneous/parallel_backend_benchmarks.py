"""
=======================================================================
Comparing parallelization with context manager using different backends
=======================================================================

In this notebook we demo how to use a context manager to change the default
backend used by the joblib implementation inside any scikit-learn object that
has a parameter `n_jobs`.

In practice, whether increasing the number of workers is helpful at improving
runtime depends on many factors. It is usually better to experiment rather than
assuming that it is always a good thing. In some cases it can be highly
detrimental to performance to run multiple copies of some estimators or
functions in parallel (see :ref:`oversubscription`).

To run this notebook locally, you may need to install dask and ray. These
packages can be installed with `pip install dask-ml "ray[default]"`.

For more information see the :ref:`User Guide <parallelism>`.

"""

# %%
import sys
import joblib
import loky
import sklearn

try:
    import dask
except ImportError:
    print("The package 'dask' is required to run this example.")
    sys.exit()

try:
    import ray
except ImportError:
    print("The package 'ray' is required to run this example.")
    sys.exit()

print(f"scikit-learn version: {sklearn.__version__}")
print(f"joblib version: {joblib.__version__}")
print(f"dask version: {dask.__version__}")
print(f"ray version: {ray.__version__}")
print(f"loky version: {loky.__version__}")

# %%
# Sample output::
#
#     scikit-learn version: 1.2.1
#     joblib version: 1.2.0
#     dask version: 2023.2.0
#     ray version: 2.2.0
#     loky version: 3.3.0
#
# This script automatically adapts to the maximum number of physical cores on
# the host. In the case of the present example, it was originally run on a
# laptop with 4 of them.

# %%
N_CORES = joblib.cpu_count(only_physical_cores=True)
print(f"number of physical cores: {N_CORES}")

# %%
# Sample output::
#
#     number of physical cores: 4
#
# We build a classification task using
# :class:`~sklearn.datasets.make_classification` and cross-validate an
# :class:`~sklearn.ensemble.HistGradientBoostingClassifier` with default
# parameters on top of it.

# %%
from sklearn.datasets import make_classification
from sklearn.ensemble import HistGradientBoostingClassifier

X, y = make_classification(n_samples=1_000, random_state=0)
hist_gbdt = HistGradientBoostingClassifier(random_state=0)

# %%
# A common setting is to estimate the performance of the model via
# cross-validation. It is sometimes interesting to set a high number of splits
# `n_splits` to improve a model's analysis. As a consequence, more computional
# resourses are needed.

# %%
from sklearn.model_selection import ShuffleSplit

cv = ShuffleSplit(n_splits=100, random_state=0)

# %%
# The computional time can still be reduced by optimizing the number of CPUs
# used via the parameter `n_jobs`.
#
# In the case of the `~sklearn.model_selection.cross_validate` function,
# `n_jobs=None` sets the number of workers in a :obj:`joblib.parallel_backend`
# context. For such function, the parallelization consists in training the
# estimator and computing the score in parallel over the cross-validation
# splits.

# %%
from time import time
from joblib import parallel_backend
from sklearn.model_selection import cross_validate


def bench(jobs, backend):
    durations = []
    msg = f"Benchmarking on {backend}:"
    print(f"\n{msg}\n" + str("-" * len(msg)))
    for n_jobs in jobs:
        with parallel_backend(backend, n_jobs=n_jobs):
            tic = time()
            cross_validate(hist_gbdt, X, y, cv=cv)
            toc = time()
        durations.append(toc - tic)
        print(f"n_jobs: {n_jobs:<3} elapsed time: {toc - tic:.3f} sec")
    return durations


# %%
from ray.util.joblib import register_ray
from dask.distributed import Client

client = Client(processes=False)  # init local dask client
ray.shutdown()  # in case there is a previously open ray session
ray.init(num_cpus=32)
register_ray()

# %%
# The scikit-learn parallelization API relies on the `loky` backend, as it is
# joblib's default backend. Here we additionally benchmark on the `threading`,
# `dask` and `ray` backends.

# %%
import numpy as np

jobs = [1, 2, 3, 4, 5, 8, 16, 32]
results = []

for backend in ["loky", "threading", "dask", "ray"]:
    durations = bench(jobs, backend)
    results.append(
        dict(
            backend=backend,
            durations=np.array(durations),
        )
    )

# %%
# Sample output::
#
#     Benchmarking on loky:
#     ---------------------
#     n_jobs: 1   elapsed time: 45.015 sec
#     n_jobs: 2   elapsed time: 17.751 sec
#     n_jobs: 3   elapsed time: 11.140 sec
#     n_jobs: 4   elapsed time: 9.403 sec
#     n_jobs: 5   elapsed time: 9.384 sec
#     n_jobs: 8   elapsed time: 10.896 sec
#     n_jobs: 16  elapsed time: 8.808 sec
#     n_jobs: 32  elapsed time: 11.410 sec
#
#     Benchmarking on threading:
#     --------------------------
#     n_jobs: 1   elapsed time: 45.032 sec
#     n_jobs: 2   elapsed time: 106.694 sec
#     n_jobs: 3   elapsed time: 98.651 sec
#     n_jobs: 4   elapsed time: 97.991 sec
#     n_jobs: 5   elapsed time: 107.621 sec
#     n_jobs: 8   elapsed time: 111.208 sec
#     n_jobs: 16  elapsed time: 114.035 sec
#     n_jobs: 32  elapsed time: 112.787 sec
#
#     Benchmarking on dask:
#     ---------------------
#     n_jobs: 1   elapsed time: 40.938 sec
#     n_jobs: 2   elapsed time: 19.351 sec
#     n_jobs: 3   elapsed time: 15.012 sec
#     n_jobs: 4   elapsed time: 12.925 sec
#     n_jobs: 5   elapsed time: 12.968 sec
#     n_jobs: 8   elapsed time: 11.780 sec
#     n_jobs: 16  elapsed time: 14.506 sec
#     n_jobs: 32  elapsed time: 20.232 sec
#
#     Benchmarking on ray:
#     --------------------
#     n_jobs: 1   elapsed time: 110.725 sec
#     n_jobs: 2   elapsed time: 111.477 sec
#     n_jobs: 3   elapsed time: 110.871 sec
#     n_jobs: 4   elapsed time: 110.623 sec
#     n_jobs: 5   elapsed time: 110.852 sec
#     n_jobs: 8   elapsed time: 3133.999 sec
#     n_jobs: 16  elapsed time: 89.752 sec
#     n_jobs: 32  elapsed time: 95.225 sec

# %%
# One can additionally plot the speedup as a function of the number of workers
# for each backend. This is not shown in the present example to avoid excessive
# computing times.

# %%
import seaborn as sns

from matplotlib import pyplot as plt

ax = sns.lineplot(x=jobs, y=jobs, color="black", label="linear growth")

for result in results:
    backend, durations = list(result.values())
    speedup = durations[0] / durations
    label = f"{backend}"
    sns.lineplot(x=jobs, y=speedup, marker="o", ax=ax, label=label)

ax.set(xscale="log", yscale="log")
ax.set_xticks(jobs)
ax.set_xticklabels(jobs)
ax.set_xlabel("number of jobs")
ax.set_ylabel("speedup")
ax.set_title("Speedup by backend and task type")
_ = plt.legend()
