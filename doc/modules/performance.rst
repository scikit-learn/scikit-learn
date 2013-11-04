.. _performance:

===========
Performance
===========

For some applications the performance (mainly speed and throughput) of
estimators is crucial. We will review here the orders of magnitude you can
expect from a number of scikit-learn estimators in different contexts and
provide some tips and tricks for overcoming performance bottlenecks.

Prediction Speed
================

One of the most straight-forward concerns one may have when using/choosing a
machine learning toolkit is the speed at which predictions can be made in a
production environment.

The main factors that influence the prediction speed are
  1. Number of features
  2. Input data representation and sparsity
  3. Model complexity
  4. Feature extraction

A last major parameter is also the possibility to do predictions in bulk or
one-at-a-time mode.

Bulk versus Atomic mode
-----------------------

In general doing predictions in bulk (many instances at the same time) is
more efficient for a number of reasons (branching predictability, CPU cache,
linear algebra libraries optimizations etc.). Here we see on a setting
with few features that independently of estimator choice the bulk mode is
always faster by 2 orders of magnitude:

.. |atomic_prediction_speed| image::  ../auto_examples/applications/images/plot_prediction_latency_1.png
    :target: ../auto_examples/applications/plot_prediction_latency.html
    :scale: 80

.. centered:: |atomic_prediction_speed|

.. |bulk_prediction_speed| image::  ../auto_examples/applications/images/plot_prediction_latency_2.png
    :target: ../auto_examples/applications/plot_prediction_latency.html
    :scale: 80

.. centered:: |bulk_prediction_speed|

To benchmark different estimators for your case you can simply change the
`n_features` parameter according to your case in this example:
:ref:`example_applications_plot_prediction_latency.py`. This should give you
an estimate of the order of magnitude of the prediction speed for your case.

Influence of the number of Features
-----------------------------------

Obviously when the number of features increases so does the memory
consumption of each example. Indeed, for a matrix of `M` instances with `N`
features, the space complexity is in `O(N.M)`. From a computing perspective
it also means that the number of basic operations (e.g. multiplications for
vector-matrix products in linear models) increases too. Here is a graph of
the evolution of the prediction speed with the number of features:

.. |influence_of_n_features_on_speed| image::  ../auto_examples/applications/images/plot_prediction_latency_3.png
    :target: ../auto_examples/applications/plot_prediction_latency.html
    :scale: 80

.. centered:: |influence_of_n_features_on_speed|

Overall you can expect the prediction time to increase at least linearly with
the number of features (non-linear cases can happen depending on the global
memory footprint and estimator).

Influence of the Input Data Representation
------------------------------------------

Numpy / Scipy support sparse matrix formats which are optimized for storing
sparse data. The main feature of sparse formats is that you don't store zeros
so if your data is sparse then you use much less memory. A non-zero value in
a sparse (`CSR or CSC <http://docs.scipy.org/doc/scipy/reference/sparse.html>`_)
representation will only take on average one 32bit integer position + the 64
bit floating point value. Using sparse input on a dense (or sparse) linear
model can speedup prediction prediction by quite a bit as only the non zero
valued features impact the dot product and thus the model predictions. Hence
if you have 100 non zeros in 1e6 dimensional space, you only need 100 multiply
+ add operation instead of 1e6.

Note that dense / dense operations benefit from both BLAS-provided SSE
vectorized operations and multithreading and lower CPU cache miss rates. Sparse
dot product is more hit or miss and does not leverage the optimized BLAS
benefit. So the sparsity should typically be quite high (10% non-zeros max,
to be checked depending on the hardware) for the sparse input representation
to be faster that the dense input representation on a machine with many CPU and
an optimized BLAS implementation.

Here is a sample code to test the sparsity of your input:

    >>> import numpy as np
    >>> def sparsity_ratio(X):
    >>>     return np.count_nonzero(X) / float(X.shape[0] * X.shape[1])
    >>> print("input sparsity ratio:", sparsity_ratio(X))

Now if you want to try to leverage sparsity for your input data you should
either build your input matrix in the CSR or CSC or call the ``to_csr()``
method or the ``csr_matrix()`` helper function from Scipy.

Prediction Throughput
=====================

Another important metric to care about when sizing production systems is the
throughput i.e. the number of predictions you can make in a given amount of
time. Here is a benchmark from the
:ref:`example_applications_plot_prediction_latency.py` example that measures
this quantity for a number of estimators on synthetic data:

.. |throughput_benchmark| image::  ../auto_examples/applications/images/plot_prediction_latency_4.png
    :target: ../auto_examples/applications/plot_prediction_latency.html
    :scale: 80

.. centered:: |throughput_benchmark|

These throughputs are achieved on a single process. An obvious way to
increase the throughput of your application is to spawn additional instances
(usually processes in Python because of the
`GIL <https://wiki.python.org/moin/GlobalInterpreterLock>`_) that share the
same model. One might also add machines to spread the load. A detailed
explanation on how to achieve this is beyond the scope of this documentation
though.

Feature Extraction Speed
========================

In many real world applications the feature extraction process (i.e. turning
raw data like database rows or network packets into numpy arrays) governs the
overall prediction time. For example here on the Reuters text classification
task the vectorization that includes parsing SGML files, tokenizing the text
and hashing it into a common vector space is taking 5 to 30 times more time
than the actual prediction code, depending on the chosen model.

 .. |computation_time| image::  ../auto_examples/applications/images/plot_out_of_core_classification_3.png
    :target: ../auto_examples/applications/plot_out_of_core_classification.html
    :scale: 80

.. centered:: |computation_time|

In many cases it is thus recommended to carefully time and profile your
feature extraction code as it may be a good place to start optimizing when
your overall speed is too slow for your application. If needed,
you can consider rewriting the feature extraction part in a lower-level,
compiled language to further speed up the overall process. The fact that
most scikit-learn models are implemented using Cython and optimized,
compiled computing libraries under the hood make them usually pretty fast.
So optimizing the feature extraction step while keeping the prediction in
python with scikit-learn estimators is usually a good way to go as it allows
for easy experimentation on the modeling side without sacrificing performance.

Tips and Tricks
===============

Linear algebra libraries
------------------------

As scikit-learn relies heavily on Numpy/Scipy and linear algebra in general it
makes sense to take explicit care of the versions of these libraries.
Basically, you ought to make sure that Numpy is built using an optimized `BLAS
<http://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms>`_ /
`LAPACK <http://en.wikipedia.org/wiki/LAPACK>`_ library.

Not all models benefit from optimized BLAS and Lapack implementations. For
instance models based on (randomized) decision trees typically do not rely on
BLAS calls in their inner loops. So do models implemented in third party C++
library (like `LinearSVC`, `LogisticRegression` from `liblinear` and SVC / SVR
from `libsvm`). On the other hand linear model implemented with a BLAS DGEMM
call (via numpy.dot) will typically benefit hugely from a tuned BLAS
implementation and lead to orders of magnitude speedup over a non-optimized
BLAS.

You can display the BLAS / LAPACK implementation used by your NumPy / SciPy /
scikit-learn install with the following commands:


    >>> from numpy.distutils.system_info import get_info
    >>> print(get_info('blas_opt'))
    >>> print(get_info('lapack_opt'))


Optimized BLAS / LAPACK implementations include:
 - Atlas (need hardware specific tuning by rebuilding on the target machine)
 - OpenBLAS
 - MKL
 - Apple Accelerate and vecLib frameworks (OSX only)

More information can be found on the `Scipy install page <http://docs.scipy
.org/doc/numpy/user/install.html>`_
and in this
`blog post <http://danielnouri.org/notes/2012/12/19/libblas-and-liblapack-issues-and-speed,-with-scipy-and-ubuntu/>`_
from Daniel Nouri which has some nice step by step install instructions.

Model Compression
-----------------

Model compression in scikit-learn only concerns linear models for the moment.
In this context it means that we want to control the model sparsity (i.e. the
number of non-zero coordinates in the model vectors). It is generally a good
idea to combine model sparsity with sparse input data representation.

Here is a sample code that illustrates the use of the ``sparsify()`` method:

    >>> clf = SGDRegressor(penalty='l1')
    >>> clf.fit(X_train, y_train)
    >>> clf.sparsify()
    >>> clf.predict(X_test)

A typical benchmark (:ref:`benchmarks_bench_sparsify.py`) on synthetic data
yields a >30% decrease in latency when both the model and input are sparsed
(with 0.000024 and 0.027400 non-zero coefficients ratio respectively).
Your mileage may vary depending on the sparsity and size of your data and
model.

Links
-----

  - `scikit-learn developer performance documentation <http://scikit-learn.org/stable/developers/performance.html>`_
  - `Scipy sparse matrix formats documentation <http://docs.scipy.org/doc/scipy/reference/sparse.html>`_
