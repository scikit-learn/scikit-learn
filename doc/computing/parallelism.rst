Parallelism and resource management
===================================

.. _parallelism:

Parallelism
-----------

Some scikit-learn estimators and utilities parallelize costly operations
using multiple CPU cores.

Depending on the type of estimator and sometimes the values of the
constructor parameters, this is either done:

- with higher-level parallelism via `joblib <https://joblib.readthedocs.io/en/latest/>`_.
- with lower-level parallelism via OpenMP, used in C or Cython code.
- with lower-level parallelism via BLAS, used by NumPy and SciPy for generic operations
  on arrays.

The ``n_jobs`` parameters of estimators always controls the amount of parallelism
managed by joblib (processes or threads depending on the joblib backend).
The thread-level parallelism managed by OpenMP in scikit-learn's own Cython code
or by BLAS & LAPACK libraries used by NumPy and SciPy operations used in scikit-learn
is always controlled by environment variables or ``threadpoolctl`` as explained below.
Note that some estimators can leverage all three kinds of parallelism at different
points of their training and prediction methods.

We describe these 3 types of parallelism in the following subsections in more detail.

Higher-level parallelism with joblib
....................................

When the underlying implementation uses joblib, the number of workers (threads
or processes) that are spawned in parallel can be controlled via the ``n_jobs``
parameter. If you choose to enable parallelism, a good starting value for
``n_jobs`` is ``-1``, which means "use all available cores." In joblib 1.7 and
later this will take into account nested thread pools, so available cores will
be smaller in worker pools.

Joblib is able to support both multi-processing and multi-threading. Whether
joblib chooses to spawn a thread or a process depends on the **backend**
that it's using.

scikit-learn generally relies on the ``loky`` backend, which is joblib's
default backend. Loky is a multi-processing backend. When doing
multi-processing, in order to avoid duplicating the memory in each process
(which isn't reasonable with big datasets), joblib will create a `memmap
<https://docs.scipy.org/doc/numpy/reference/generated/numpy.memmap.html>`_
that all processes can share, when the data is bigger than 1MB.

In some specific cases (when the code that is run in parallel releases the
GIL), scikit-learn will indicate to ``joblib`` that a multi-threading
backend is preferable. On free-threaded Python, which doesn't have a GIL,
the threading backend will be enabled in additional cases.

On free-threaded Python (sometimes identified by a "t" suffix on the version,
e.g. "3.15t"), the threaded backend will often be faster than ``loky``.

As a user, you may control the backend that joblib will use (regardless of
what scikit-learn recommends) by using a context manager::

    from joblib import parallel_backend

    with parallel_backend('threading', n_jobs=-1):
        # Your scikit-learn code here

Please refer to the `joblib's docs
<https://joblib.readthedocs.io/en/latest/parallel.html#thread-based-parallelism-vs-process-based-parallelism>`_
for more details.

In practice, whether parallelism is helpful at improving runtime depends on
many factors. It is usually a good idea to experiment rather than assuming
that increasing the number of workers is always a good thing. In some cases
it can be highly detrimental to performance to run multiple copies of some
estimators or functions in parallel (see :ref:`oversubscription<oversubscription>` below).

.. note::

    **Startup Overhead**

    When using ``n_jobs > 1`` (or ``n_jobs=-1``) with a process pool, you may
    observe a delay the first time a parallel function is called. This is
    expected behavior caused by the overhead of starting the Python worker
    processes. Subsequent calls will be faster as they reuse the existing pool
    of workers.

.. note::

    Where (and how) parallelization happens in the estimators using joblib by
    specifying `n_jobs` is currently poorly documented.
    Please help us by improving our docs and tackle `issue 14228
    <https://github.com/scikit-learn/scikit-learn/issues/14228>`_!

.. _lower-level-parallelism-with-openmp:

Lower-level parallelism with OpenMP
...................................

OpenMP is used to parallelize code written in Cython or C, relying on
multi-threading exclusively. By default, the implementations using OpenMP
will use as many threads as possible, i.e. as many threads as logical cores.

Within a single Python thread, you can control the exact number of threads that
are used via `threadpoolctl` as explained by `this piece of documentation
<https://github.com/joblib/threadpoolctl/#setting-the-maximum-size-of-thread-pools>`_.

Within the context of thread or process pools created by ``joblib``, ``joblib``
1.7 and later will manage OpenMP limits appropriately (earlier versions only
managed process pools).

:ref:`SKLEARN_PAIRWISE_DIST_CHUNK_SIZE <envvar_SKLEARN_PAIRWISE_DIST_CHUNK_SIZE>`
tunes chunk size for OpenMP-parallel pairwise distance kernels in Cython (e.g.
nearest neighbors). This is an advanced setting.

Parallel NumPy and SciPy routines from numerical libraries
..........................................................

scikit-learn relies heavily on NumPy and SciPy, which internally call
multi-threaded linear algebra routines (BLAS & LAPACK) implemented in libraries
such as MKL, OpenBLAS or BLIS.

Starting with ``joblib`` 1.7, scikit-learn relies on ``joblib`` to set the
number of threads in these libraries, and they should be set to reasonable
values automatically.

.. _oversubscription:

Oversubscription: spawning too many threads
...........................................

It is generally recommended to avoid using significantly more processes or
threads than the number of CPUs on a machine. Over-subscription happens when a
program is running too many threads at the same time. This can happen through
two mechanisms:

* Nested algorithms that use ``joblib.Parallel`` in both the inner and outer algorithm.
* An outer ``joblib.Parallel`` whose workers then call into BLAS or OpenMP.

Starting with ``joblib`` v1.7 both are much less of a concern, though some
configuration may be required.

For nested parallel algorithms, as long as you configure the inner algorithms
with ``n_jobs=-1``, meaning "use all available cores", inner algorithms
parallelism will take into account the number of threads that the outer
algorithm uses. For example if the outer algorithm uses ``n_jobs=2`` and the
inner algorithm uses ``n_jobs=-1``, the inner algorithm will use ``number of
cores // 2`` workers. Avoid hardcoding a specific number of workers in the inner
algorithm, e.g. 4, if you're using ``n_jobs=-1`` in the outer algorithm.

For BLAS and OpenMP, scikit-learn together with ``joblib`` 1.7 or later should
automatically handle this for you.
Suppose you have a machine with 8 CPUs. Consider a case where you're running a
:class:`~sklearn.model_selection.GridSearchCV` (parallelized with joblib) with
``n_jobs=8`` over a :class:`~sklearn.ensemble.HistGradientBoostingClassifier`
(parallelized with OpenMP). Naively, you'd expect each instance of
:class:`~sklearn.ensemble.HistGradientBoostingClassifier` will spawn 8 threads
(since you have 8 CPUs). That's a total of ``8 * 8 = 64`` threads, which leads
to oversubscription of threads for physical CPU resources and thus to scheduling
overhead.

To prevent this situation, starting from ``joblib >= 1.7``, when the ``loky``
backend (the default) or the ``threading`` backend is used, joblib will tell its
child processes and threads to limit the number of threads they can use, so as
to avoid oversubscription. (In older versions only ``loky`` is restricted, so
you will get oversubscription when using threads, so upgrading is recommended if
you're on an older version.) In practice the heuristic that joblib uses is to
tell the processes to use ``max_threads = n_cpus // n_jobs``. Back to our
example from above, since the joblib backend of
:class:`~sklearn.model_selection.GridSearchCV` is ``loky``, each process will
only be able to use 1 thread instead of 8, thus mitigating the oversubscription
issue.

Note that:

- All scikit-learn estimators that explicitly rely on OpenMP in their Cython code
  always use ``threadpoolctl`` internally to automatically adapt the numbers of
  threads used by OpenMP and potentially nested BLAS calls so as to avoid
  oversubscription.

You will find additional details about joblib mitigation of oversubscription
in `joblib documentation
<https://joblib.readthedocs.io/en/latest/parallel.html#avoiding-over-subscription-of-cpu-resources>`_.

You will find additional details about parallelism in numerical python libraries
in `this document from Thomas J. Fan
<https://thomasjpfan.github.io/parallelism-python-libraries-design/>`_; some of
the details may be out of date, given subsequent work on ``joblib`` and
``threadpoolctl``.
