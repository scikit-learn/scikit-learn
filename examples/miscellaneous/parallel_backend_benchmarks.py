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

You may need to install dask and ray to run this notebook. These packages can be
installed with `pip install dask "ray[default]"`.

For more information on parallelism, see the :ref:`User Guide <parallelism>`.

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

# %%
# This script was originally run using the following versions for python and the
# relevant packages:

# %%
print(f"python version: {sys.version.split(' ')[0]}")
print(f"scikit-learn version: {sklearn.__version__}")
print(f"joblib version: {joblib.__version__}")
print(f"dask version: {dask.__version__}")
print(f"ray version: {ray.__version__}")
print(f"loky version: {loky.__version__}")

# %%
# Sample output::
#
#     python version: 3.9.16
#     scikit-learn version: 1.2.1
#     joblib version: 1.2.0
#     dask version: 2023.2.0
#     ray version: 2.2.0
#     loky version: 3.3.0
#
# This script also automatically adapts to the maximum number of physical cores
# on the host. In the case of the present example, it was originally run on a
# laptop with 4 of them.

# %%
N_CORES = joblib.cpu_count(only_physical_cores=True)
print(f"number of physical cores: {N_CORES}")

# %%
# Sample output::
#
#     number of physical cores: 4
#
# Once settled the specification, we build a classification task using
# :class:`~sklearn.datasets.make_classification` and cross-validate an
# :class:`~sklearn.ensemble.HistGradientBoostingClassifier` with default
# parameters on top of it.

# %%
from sklearn.datasets import make_classification
from sklearn.ensemble import HistGradientBoostingClassifier

X, y = make_classification(n_samples=1_000, random_state=0)
clf = HistGradientBoostingClassifier(random_state=0)

# %%
# A common setting is to estimate the performance of the model via
# cross-validation. It is sometimes interesting to set a high number of splits
# `n_splits` to improve a model's analysis. As a consequence, more computional
# resourses are needed.

# %%
from sklearn.model_selection import ShuffleSplit

cv = ShuffleSplit(n_splits=10, random_state=0)

# %%
# The computional time can still be reduced by optimizing the number of CPUs
# used via the parameter `n_jobs`.
#
# In the case of the `~sklearn.model_selection.cross_validate` function, the
# default `n_jobs=None` allows us to set the number of workers within a
# :obj:`joblib.parallel_backend` context manager. For such function, the
# parallelization consists in training the estimator and computing the score in
# parallel over the cross-validation splits.

# %%
from time import time
from sklearn.model_selection import cross_validate
from threadpoolctl import threadpool_limits
import numpy as np

n_threads_grid = 2 ** np.arange(np.log2(2 * N_CORES).astype(np.int32) + 1)

for n_threads in n_threads_grid:
    tic = time()
    with threadpool_limits(limits=int(n_threads)):
        cross_validate(clf, X, y, cv=cv)
    toc = time()
    print(f"n_threads: {n_threads}, elapsed time: {toc - tic:.3f} sec")

# %%
from joblib import parallel_backend


def bench(n_jobs_grid, backend):
    durations = []
    msg = f"Benchmarking on {backend}:"
    print(f"\n{msg}\n" + str("-" * len(msg)))
    for n_jobs in n_jobs_grid:
        with parallel_backend(backend, n_jobs=int(n_jobs)):
            tic = time()
            cross_validate(clf, X, y, cv=cv)
            toc = time()
        durations.append(toc - tic)
        print(f"n_jobs: {n_jobs:<3} elapsed time: {toc - tic:.3f} sec")
    return durations


# %%
# The scikit-learn parallelization API relies on the `loky` backend, as it is
# joblib's default backend. Here we additionally benchmark on the `threading`,
# `dask` and `ray` backends. The last two require to be init as follows:

# %%
from ray.util.joblib import register_ray
from dask.distributed import Client

client = Client(processes=False)  # init dask client
ray.shutdown()  # in case there is a previously open ray session
ray.init(num_cpus=N_CORES)  # init ray client
register_ray()

# %%
# We define a grid of the number of workers spaced in powers of 2. To avoid
# oversubscription, the grid's maximal value is set to be `N_CORES`.

# %%
n_jobs_grid = 2 ** np.arange(np.log2(N_CORES).astype(np.int32) + 1)
results = []

for backend in ["loky", "threading", "dask", "ray"]:
    durations = bench(n_jobs_grid, backend)
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
#     n_jobs: 4   elapsed time: 9.403 sec
#
#     Benchmarking on threading:
#     --------------------------
#     n_jobs: 1   elapsed time: 45.032 sec
#     n_jobs: 2   elapsed time: 106.694 sec
#     n_jobs: 4   elapsed time: 97.991 sec
#
#     Benchmarking on dask:
#     ---------------------
#     n_jobs: 1   elapsed time: 40.938 sec
#     n_jobs: 2   elapsed time: 19.351 sec
#     n_jobs: 4   elapsed time: 12.925 sec
#
#     Benchmarking on ray:
#     --------------------
#     n_jobs: 1   elapsed time: 41.569 sec
#     n_jobs: 2   elapsed time: 15.271 sec
#     n_jobs: 4   elapsed time: 11.742 sec

# %%
# One can additionally plot the speedup as a function of the number of workers
# for each backend. This is not shown in the present example to avoid excessive
# computing times.

# %%
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

ax = sns.lineplot(x=n_jobs_grid, y=n_jobs_grid, color="black", label="linear growth")

for result in results:
    backend, durations = list(result.values())
    speedup = durations[0] / durations
    label = f"{backend}"
    sns.lineplot(x=n_jobs_grid, y=speedup, marker="o", ax=ax, label=label)

ax.set(xscale="log", yscale="log")
ax.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
ax.set_xticks(n_jobs_grid)
ax.set_xticklabels(n_jobs_grid)
ax.set_xlabel("number of jobs")
ax.set_ylabel("speedup")
ax.set_title("Speedup by backend and task type")
_ = plt.legend()
