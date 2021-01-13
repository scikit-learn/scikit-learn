.. Places parent toc into the sidebar

:parenttoc: True

Parallelism, resource management, and configuration
===================================================

.. _parallelism:

Parallelism
-----------

Some scikit-learn estimators and utilities can parallelize costly operations
using multiple CPU cores, thanks to the following components:

- via the `joblib <https://joblib.readthedocs.io/en/latest/>`_ library. In
  this case the number of threads or processes can be controlled with the
  ``n_jobs`` parameter.
- via OpenMP, used in C or Cython code.

In addition, some of the numpy routines that are used internally by
scikit-learn may also be parallelized if numpy is installed with specific
numerical libraries such as MKL, OpenBLAS, or BLIS.

We describe these 3 scenarios in the following subsections.

Joblib-based parallelism
........................

When the underlying implementation uses joblib, the number of workers
(threads or processes) that are spawned in parallel can be controlled via the
``n_jobs`` parameter.

.. note::

    Where (and how) parallelization happens in the estimators is currently
    poorly documented. Please help us by improving our docs and tackle `issue
    14228 <https://github.com/scikit-learn/scikit-learn/issues/14228>`_!

Joblib is able to support both multi-processing and multi-threading. Whether
joblib chooses to spawn a thread or a process depends on the **backend**
that it's using.

Scikit-learn generally relies on the ``loky`` backend, which is joblib's
default backend. Loky is a multi-processing backend. When doing
multi-processing, in order to avoid duplicating the memory in each process
(which isn't reasonable with big datasets), joblib will create a `memmap
<https://docs.scipy.org/doc/numpy/reference/generated/numpy.memmap.html>`_
that all processes can share, when the data is bigger than 1MB.

In some specific cases (when the code that is run in parallel releases the
GIL), scikit-learn will indicate to ``joblib`` that a multi-threading
backend is preferable.

As a user, you may control the backend that joblib will use (regardless of
what scikit-learn recommends) by using a context manager::

    from joblib import parallel_backend

    with parallel_backend('threading', n_jobs=2):
        # Your scikit-learn code here

Please refer to the `joblib's docs
<https://joblib.readthedocs.io/en/latest/parallel.html#thread-based-parallelism-vs-process-based-parallelism>`_
for more details.

In practice, whether parallelism is helpful at improving runtime depends on
many factors. It is usually a good idea to experiment rather than assuming
that increasing the number of workers is always a good thing. In some cases
it can be highly detrimental to performance to run multiple copies of some
estimators or functions in parallel (see oversubscription below).

OpenMP-based parallelism
........................

OpenMP is used to parallelize code written in Cython or C, relying on
multi-threading exclusively. By default (and unless joblib is trying to
avoid oversubscription), the implementation will use as many threads as
possible.

You can control the exact number of threads that are used via the
``OMP_NUM_THREADS`` environment variable:

.. prompt:: bash $

    OMP_NUM_THREADS=4 python my_script.py

Parallel Numpy routines from numerical libraries
................................................

Scikit-learn relies heavily on NumPy and SciPy, which internally call
multi-threaded linear algebra routines implemented in libraries such as MKL,
OpenBLAS or BLIS.

The number of threads used by the OpenBLAS, MKL or BLIS libraries can be set
via the ``MKL_NUM_THREADS``, ``OPENBLAS_NUM_THREADS``, and
``BLIS_NUM_THREADS`` environment variables.

Please note that scikit-learn has no direct control over these
implementations. Scikit-learn solely relies on Numpy and Scipy.

.. note::
    At the time of writing (2019), NumPy and SciPy packages distributed on
    pypi.org (used by ``pip``) and on the conda-forge channel are linked
    with OpenBLAS, while conda packages shipped on the "defaults" channel
    from anaconda.org are linked by default with MKL.


Oversubscription: spawning too many threads
...........................................

It is generally recommended to avoid using significantly more processes or
threads than the number of CPUs on a machine. Over-subscription happens when
a program is running too many threads at the same time.

Suppose you have a machine with 8 CPUs. Consider a case where you're running
a :class:`~sklearn.model_selection.GridSearchCV` (parallelized with joblib)
with ``n_jobs=8`` over a
:class:`~sklearn.ensemble.HistGradientBoostingClassifier` (parallelized with
OpenMP). Each instance of
:class:`~sklearn.ensemble.HistGradientBoostingClassifier` will spawn 8 threads
(since you have 8 CPUs). That's a total of ``8 * 8 = 64`` threads, which
leads to oversubscription of physical CPU resources and to scheduling
overhead.

Oversubscription can arise in the exact same fashion with parallelized
routines from MKL, OpenBLAS or BLIS that are nested in joblib calls.

Starting from ``joblib >= 0.14``, when the ``loky`` backend is used (which
is the default), joblib will tell its child **processes** to limit the
number of threads they can use, so as to avoid oversubscription. In practice
the heuristic that joblib uses is to tell the processes to use ``max_threads
= n_cpus // n_jobs``, via their corresponding environment variable. Back to
our example from above, since the joblib backend of
:class:`~sklearn.model_selection.GridSearchCV` is ``loky``, each process will
only be able to use 1 thread instead of 8, thus mitigating the
oversubscription issue.

Note that:

- Manually setting one of the environment variables (``OMP_NUM_THREADS``,
  ``MKL_NUM_THREADS``, ``OPENBLAS_NUM_THREADS``, or ``BLIS_NUM_THREADS``)
  will take precedence over what joblib tries to do. The total number of
  threads will be ``n_jobs * <LIB>_NUM_THREADS``. Note that setting this
  limit will also impact your computations in the main process, which will
  only use ``<LIB>_NUM_THREADS``. Joblib exposes a context manager for
  finer control over the number of threads in its workers (see joblib docs
  linked below).
- Joblib is currently unable to avoid oversubscription in a
  multi-threading context. It can only do so with the ``loky`` backend
  (which spawns processes).

You will find additional details about joblib mitigation of oversubscription
in `joblib documentation
<https://joblib.readthedocs.io/en/latest/parallel.html#avoiding-over-subscription-of-cpu-ressources>`_.


Configuration switches
-----------------------

Python runtime
..............

:func:`sklearn.set_config` controls the following behaviors:

:assume_finite:

    used to skip validation, which enables faster computations but may
    lead to segmentation faults if the data contains NaNs.

:working_memory:

    the optimal size of temporary arrays used by some algorithms.

.. _environment_variable:

Environment variables
......................

These environment variables should be set before importing scikit-learn.

:SKLEARN_SITE_JOBLIB:

    When this environment variable is set to a non zero value,
    scikit-learn uses the site joblib rather than its vendored version.
    Consequently, joblib must be installed for scikit-learn to run.
    Note that using the site joblib is at your own risks: the versions of
    scikit-learn and joblib need to be compatible. Currently, joblib 0.11+
    is supported. In addition, dumps from joblib.Memory might be incompatible,
    and you might loose some caches and have to redownload some datasets.

    .. deprecated:: 0.21

       As of version 0.21 this parameter has no effect, vendored joblib was
       removed and site joblib is always used.

:SKLEARN_ASSUME_FINITE:

    Sets the default value for the `assume_finite` argument of
    :func:`sklearn.set_config`.

:SKLEARN_WORKING_MEMORY:

    Sets the default value for the `working_memory` argument of
    :func:`sklearn.set_config`.

:SKLEARN_SEED:

    Sets the seed of the global random generator when running the tests,
    for reproducibility.

:SKLEARN_SKIP_NETWORK_TESTS:

    When this environment variable is set to a non zero value, the tests
    that need network access are skipped. When this environment variable is
    not set then network tests are skipped.
