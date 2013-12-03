.. _computational_performance:

=========================
Computational Performance
=========================

For some applications the performance (mainly latency and throughput at
prediction time) of estimators is crucial. It may also be of interest to
consider the training throughput but this is often less important in a
production setup (where it often takes place offline).

We will review here the orders of magnitude you can expect from a number of
scikit-learn estimators in different contexts and provide some tips and
tricks for overcoming performance bottlenecks.

Prediction latency is measured as the elapsed time necessary to make a
prediction (e.g. in micro-seconds). Latency is often viewed as a distribution
and operations engineers often focus on the latency at a given percentile of
this distribution (e.g. the 90 percentile).

Prediction throughput is defined as the number of predictions the software can
deliver in a given amount of time (e.g. in predictions per second).

An important aspect of performance optimization is also that it can hurt
prediction accuracy. Indeed, simpler models (e.g. linear instead of
non-linear, or with fewer parameters) often run faster but are not always able
to take into account the same exact properties of the data as more complex ones.

Prediction Latency
==================

One of the most straight-forward concerns one may have when using/choosing a
machine learning toolkit is the latency at which predictions can be made in a
production environment.

The main factors that influence the prediction latency are
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
always faster, and for some of them by 1 to 2 orders of magnitude:

.. |atomic_prediction_latency| image::  ../auto_examples/applications/images/plot_prediction_latency_1.png
    :target: ../auto_examples/applications/plot_prediction_latency.html
    :scale: 80

.. centered:: |atomic_prediction_latency|

.. |bulk_prediction_latency| image::  ../auto_examples/applications/images/plot_prediction_latency_2.png
    :target: ../auto_examples/applications/plot_prediction_latency.html
    :scale: 80

.. centered:: |bulk_prediction_latency|

To benchmark different estimators for your case you can simply change the
``n_features`` parameter according to your case in this example:
:ref:`example_applications_plot_prediction_latency.py`. This should give you
an estimate of the order of magnitude of the prediction latency for your case.

Influence of the number of Features
-----------------------------------

Obviously when the number of features increases so does the memory
consumption of each example. Indeed, for a matrix of `M` instances with `N`
features, the space complexity is in `O(N.M)`. From a computing perspective
it also means that the number of basic operations (e.g. multiplications for
vector-matrix products in linear models) increases too. Here is a graph of
the evolution of the prediction latency with the number of features:

.. |influence_of_n_features_on_latency| image::  ../auto_examples/applications/images/plot_prediction_latency_3.png
    :target: ../auto_examples/applications/plot_prediction_latency.html
    :scale: 80

.. centered:: |influence_of_n_features_on_latency|

Overall you can expect the prediction time to increase at least linearly with
the number of features (non-linear cases can happen depending on the global
memory footprint and estimator).

Influence of the Input Data Representation
------------------------------------------

Scipy support sparse matrix datastructures which are optimized for storing
sparse data. The main feature of sparse formats is that you don't store zeros
so if your data is sparse then you use much less memory. A non-zero value in
a sparse (`CSR or CSC <http://docs.scipy.org/doc/scipy/reference/sparse.html>`_)
representation will only take on average one 32bit integer position + the 64
bit floating point value + an additional 32bit per row or column in the matrix.
Using sparse input on a dense (or sparse) linear model can speedup prediction
by quite a bit as only the non zero valued features impact the dot product
and thus the model predictions. Hence if you have 100 non zeros in 1e6
dimensional space, you only need 100 multiply + add operation instead of 1e6.

Note that dense / dense operations benefit from both BLAS-provided SSE
vectorized operations and multithreading and lower CPU cache miss rates. Sparse
dot product is more hit or miss and does not leverage the optimized BLAS
benefit. So the sparsity should typically be quite high (10% non-zeros max,
to be checked depending on the hardware) for the sparse input representation
to be faster than the dense input representation on a machine with many CPU and
an optimized BLAS implementation.

Here is a sample code to test the sparsity of your input:

    >>> from sklearn.utils.fixes import count_nonzero
    >>> def sparsity_ratio(X):
    >>>     return 1.0 - count_nonzero(X) / float(X.shape[0] * X.shape[1])
    >>> print("input sparsity ratio:", sparsity_ratio(X))

As a rule of thumb you can consider that if the sparsity ratio is greater
than 90% you can probably benefit from sparse formats. Now if you want to try
to leverage sparsity for your input data you should either build your input
matrix in the CSR or CSC or call the ``to_csr()`` method or the ``csr_matrix()``
helper function from Scipy.

Influence of the Model Complexity
---------------------------------

Generally speaking, when model complexity increases, predictive power and
latency are supposed to increase. Increasing predictive power is usually
interesting, but for many applications we would better not increase
prediction latency too much. We will now review this idea for different
families of supervised models.

For linear models (e.g. Lasso, ElasticNet, SGDClassifier/Regressor,
Ridge & RidgeClassifier, PassiveAgressiveClassifier/Regressor, LinearSVC,
LogisticRegression...) the decision function that is applied at prediction
time is the same, so latency should be equivalent. Of course the particular
values (and sparsity) will change depending on how the model was trained but
the type of operation is the same (a dot product).

Here is an example using
:class:`sklearn.linear_model.stochastic_gradient.SGDClassifier` with the
``elasticnet`` penalty. The regularization power is globally controlled by
the ``alpha`` parameter. With a sufficiently high ``alpha``,
one can then play with the ``l1_ratio`` parameter of ``elasticnet`` to
enforce various levels of sparsity in the model coefficients. Higher sparsity
here is interpreted as less model complexity as we need less coefficients to
describe it fully. Of course sparsity influences in turn the prediction time
as the sparse dot-product takes time roughly proportional to the number of
non-zero coefficients.

.. |en_model_complexity| image::  ../auto_examples/applications/images/plot_model_complexity_influence_1.png
    :target: ../auto_examples/applications/plot_model_complexity_influence.html
    :scale: 80

.. centered:: |en_model_complexity|

For the SVM family of algorithms with a non-linear kernel, the latency is tied
to the number of support vectors (the fewer the faster). Latency and
throughput should (asymptotically) grow linearly with the number of support
vectors in a SVC or SVR model. The kernel will also influence the latency as
it is used to compute the projection of the input vector once per support
vector. In the following graph the ``nu`` parameter of
:class:`sklearn.svm.classes.NuSVR` was used to influence the number of
support vectors.

.. |nusvr_model_complexity| image::  ../auto_examples/applications/images/plot_model_complexity_influence_2.png
    :target: ../auto_examples/applications/plot_model_complexity_influence.html
    :scale: 80

.. centered:: |nusvr_model_complexity|

For ensemble of trees (e.g. RandomForest, GBT, ExternalTrees etc) the
number of trees and their depth play the most important role. Latency and
throughput should scale linearly with the number of trees. In this case
we used directly the ``n_estimators`` parameter of
:class:`sklearn.ensemble.gradient_boosting.GradientBoostingRegressor`.

.. |gbt_model_complexity| image::  ../auto_examples/applications/images/plot_model_complexity_influence_3.png
    :target: ../auto_examples/applications/plot_model_complexity_influence.html
    :scale: 80

.. centered:: |gbt_model_complexity|

In any case be warned that playing with model complexity can hurt accuracy as
mentionned above. For instance a non-linearly separable problem can be dealt
with a speedy linear model but prediction power will very likely suffer in
the process.

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

Feature Extraction Latency
==========================

In many real world applications the feature extraction process (i.e. turning
raw data like database rows or network packets into numpy arrays) governs the
overall prediction time. For example here on the Reuters text classification
task the whole preparation that includes reading and parsing SGML files,
tokenizing the text and hashing it into a common vector space is taking 100
to 500 times more time than the actual prediction code, depending on the chosen
model.

 .. |prediction_time| image::  ../auto_examples/applications/images/plot_out_of_core_classification_4.png
    :target: ../auto_examples/applications/plot_out_of_core_classification.html
    :scale: 80

.. centered:: |prediction_time|

In many cases it is thus recommended to carefully time and profile your
feature extraction code as it may be a good place to start optimizing when
your overall latency is too slow for your application. If needed,
you can consider rewriting the feature extraction part in a lower-level,
compiled language to further speed up the overall process. Most scikit-learn
models are usually pretty fast as they are implemented either with compiled
Cython extensions or optimized computing libraries. So optimizing the feature
extraction step while keeping the prediction in python with scikit-learn
estimators is usually a good way to go as it allows for easy experimentation
on the modeling side without sacrificing performance.

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
library (like ``LinearSVC``, ``LogisticRegression`` from ``liblinear`` and SVC /
SVR from ``libsvm``). On the other hand a linear model implemented with a BLAS
DGEMM call (via ``numpy.dot``) will typically benefit hugely from a tuned BLAS
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
from Daniel Nouri which has some nice step by step install instructions for
Debian / Ubuntu.

Model Compression
-----------------

Model compression in scikit-learn only concerns linear models for the moment.
In this context it means that we want to control the model sparsity (i.e. the
number of non-zero coordinates in the model vectors). It is generally a good
idea to combine model sparsity with sparse input data representation.

Here is a sample code that illustrates the use of the ``sparsify()`` method:

    >>> clf = SGDRegressor(penalty='elasticnet', l1_ratio=0.25)
    >>> clf.fit(X_train, y_train).sparsify()
    >>> clf.predict(X_test)

In this example we prefer the ``elasticnet`` penalty as it is often a good
compromise between model compactness and prediction power. One can also
further tune the ``l1_ratio`` parameter (in combination with the
regularization strength ``alpha``) to control this tradeoff.

A typical `benchmark <https://github.com/scikit-learn/scikit-learn/tree/masternchmarks/bench_sparsify.py>`_
on synthetic data yields a >30% decrease in latency when both the model and
input are sparse (with 0.000024 and 0.027400 non-zero coefficients ratio
respectively). Your mileage may vary depending on the sparsity and size of
your data and model.
Furthermore, sparsifying can be very useful to reduce the memory usage of
predictive models deployed on production servers.

Links
-----

  - `scikit-learn developer performance documentation <../developers/performance.html>`_
  - `Scipy sparse matrix formats documentation <http://docs.scipy.org/doc/scipy/reference/sparse.html>`_
