.. Places parent toc into the sidebar

:parenttoc: True

.. _randomness:

Controlling randomness: good practices and common pitfalls
==========================================================

Some scikit-learn objects are inherently random. These are usually estimators
(e.g. :class:`~sklearn.ensemble.RandomForestClassifier`) and cross-validation
splitters (e.g. :class:`~sklearn.model_selection.KFold`). The randomness of
these objects is controlled via their `random_state` parameter, as described
in the :term:`Glossary <random_state>`. This section expands on the glossary
entry, and describes good practices and common pitfalls w.r.t. to this
subtle parameter.

.. note:: Recommendation summary

    For an optimal statistical significance of cross-validation (CV) results,
    pass `RandomState` instances when creating estimators, or leave
    `random_state` to None. Passing integers to CV splitters is usually the
    safest option, although `RandomState` instances are fine if you know what
    you are doing. For both estimators and splitters, passing an integer vs
    passing an instance (or None) leads to subtle but significant
    differences, especially for CV procedures. These differences are
    important to understand when reporting results.
    
    For reproducible results across executions, remove any use of
    `random_state=None`.

Using None or `RandomState` instances, and repeated calls to `fit` and `split`
------------------------------------------------------------------------------

The `random_state` parameter determines whether multiple calls to :term:`fit`
(for estimators) or to :term:`split` (for CV splitters) will produce the same
results, according to these rules:

- If an integer is passed, calling `fit` or `split` multiple times always
  yields the same results.
- If `None` or a `RandomState` instance is passed: `fit` and `split` will
  yield different results each time they are called, and the succession of
  calls explores all sources of entropy. `None` is the default value for all
  `random_state` parameters.

We here illustrate these rules for both estimators and CV splitters.

.. note::
    Since passing `random_state=None` is equivalent to passing the global
    `RandomState` instance from `numpy`
    (`random_state=np.random.mtrand._rand`), we will not explicitly mention
    `None` here. Everything that applies to instances also applies to using
    `None`.

Estimators
..........

Passing instances means that calling `fit` multiple times will not yield the
same results, even if the estimator is fitted on the same data and with the
same hyper-parameters::

    >>> from sklearn.linear_model import SGDClassifier
    >>> from sklearn.datasets import make_classification
    >>> import numpy as np

    >>> rng = np.random.RandomState(0)
    >>> X, y = make_classification(n_features=5, random_state=rng)
    >>> sgd = SGDClassifier(random_state=rng)

    >>> sgd.fit(X, y).coef_
    array([[ 8.85418642,  4.79084103, -3.13077794,  8.11915045, -0.56479934]])

    >>> sgd.fit(X, y).coef_
    array([[ 6.70814003,  5.25291366, -7.55212743,  5.18197458,  1.37845099]])

We can see from the snippet above that repeatedly calling `sgd.fit` has
produced different models, even if the data was the same. This is because the
Random Number Generator (RNG) of the estimator is consumed (i.e. mutated)
when `fit` is called, and this mutated RNG will be used in the subsequent
calls to `fit`. In addition, the `rng` object is shared across all objects
that use it, and as a consequence, these objects become somewhat
inter-dependent. For example, two estimators that share the same
`RandomState` instance will influence each other, as we will see later when
we discuss cloning. This point is important to keep in mind when debugging.

If we had passed an integer to the `random_state` parameter of the
:class:`~sklearn.ensemble.RandomForestClassifier`, we would have obtained the
same models, and thus the same scores each time. When we pass an integer, the
same RNG is used across all calls to `fit`. What internally happens is that
even though the RNG is consumed when `fit` is called, it is always reset to
its original state at the beginning of `fit`.

CV splitters
............

Randomized CV splitters have a similar behavior when a `RandomState`
instance is passed; calling `split` multiple times yields different data
splits::

    >>> from sklearn.model_selection import KFold
    >>> import numpy as np

    >>> X = y = np.arange(10)
    >>> rng = np.random.RandomState(0)
    >>> cv = KFold(n_splits=2, shuffle=True, random_state=rng)

    >>> for train, test in cv.split(X, y):
    ...     print(train, test)
    [0 3 5 6 7] [1 2 4 8 9]
    [1 2 4 8 9] [0 3 5 6 7]

    >>> for train, test in cv.split(X, y):
    ...     print(train, test)
    [0 4 6 7 8] [1 2 3 5 9]
    [1 2 3 5 9] [0 4 6 7 8]

We can see that the splits are different from the second time `split` is
called. This may lead to unexpected results if you compare the performance of
multiple estimators by calling `split` many times, as we will see in the next
section.

Common pitfalls and subtleties
------------------------------

While the rules that govern the `random_state` parameter are seemingly simple,
they do however have some subtle implications. In some cases, this can even
lead to wrong conclusions.

Estimators
..........

**Different `random_state` types lead to different cross-validation
procedures**

Depending on the type of the `random_state` parameter, estimators will behave
differently, especially in cross-validation procedures. Consider the
following snippet::

    >>> from sklearn.ensemble import RandomForestClassifier
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.model_selection import cross_val_score
    >>> import numpy as np

    >>> X, y = make_classification(random_state=0)

    >>> rf_123 = RandomForestClassifier(random_state=123)
    >>> cross_val_score(rf_123, X, y)
    array([0.85, 0.95, 0.95, 0.9 , 0.9 ])

    >>> rf_inst = RandomForestClassifier(random_state=np.random.RandomState(0))
    >>> cross_val_score(rf_inst, X, y)
    array([0.9 , 0.95, 0.95, 0.9 , 0.9 ])

We see that the cross-validated scores of `rf_123` and `rf_inst` are
different, as should be expected since we didn't pass the same `random_state`
parameter. However, the difference between these scores is more subtle than
it looks, and **the cross-validation procedures that were performed by**
:func:`~sklearn.model_selection.cross_val_score` **significantly differ in
each case**:

- Since `rf_123` was passed an integer, every call to `fit` uses the same RNG:
  this means that all random characteristics of the random forest estimator
  will be the same for each of the 5 folds of the CV procedure. In
  particular, the (randomly chosen) subset of features of the estimator will
  be the same across all folds.
- Since `rf_inst` was passed a `RandomState` instance, each call to `fit`
  starts from a different RNG. As a result, the random subset of features
  will be different for each folds.

While having a constant estimator RNG across folds isn't inherently wrong,
results are more significant if we allow the estimator RNG to vary for each
fold, so passing an instance may be preferable (more on this later).

.. note::
    Here, :func:`~sklearn.model_selection.cross_val_score` will use a
    non-randomized CV splitter (as is the default), so both estimators will
    be evaluated on the same splits. This section is not about variability in
    the splits. Also, whether we pass an integer or an instance to
    :func:`~sklearn.datasets.make_classification` isn't relevant for our
    illustration purpose: what matters is what we pass to the
    :class:`~sklearn.ensemble.RandomForestClassifier` estimator.

**Cloning**

Another subtle side effect of passing `RandomState` instances is how
:func:`~sklearn.clone` will work::

    >>> from sklearn import clone
    >>> from sklearn.ensemble import RandomForestClassifier
    >>> import numpy as np

    >>> rng = np.random.RandomState(0)
    >>> a = RandomForestClassifier(random_state=rng)
    >>> b = clone(a)

Since a `RandomState` instance was passed to `a`, `a` and `b` are not clones
in the strict sense, but rather clones in the statistical sense: `a` and `b`
will still be different models, even when calling `fit(X, y)` on the same
data. Moreover, `a` and `b` will influence each-other since they share the
same internal RNG: calling `a.fit` will consume `b`'s RNG, and calling
`b.fit` will consume `a`'s RNG, since they are the same. This bit is true for
any estimators that share a `random_state` parameter; it is not specific to
clones.

If an integer were passed, `a` and `b` would be exact clones and they would not
influence each other.

.. warning::
    Even though :func:`~sklearn.clone` is rarely used in user code, it is
    called pervasively throughout scikit-learn codebase: in particular, most
    meta-estimators that accept non-fitted estimators call
    :func:`~sklearn.clone` internally
    (:class:`~sklearn.model_selection.GridSearchCV`,
    :class:`~sklearn.ensemble.StackingClassifier`,
    :class:`~sklearn.calibration.CalibratedClassifierCV`, etc.).

CV splitters
............

When passed a `RandomState` instance, CV splitters yield different splits
each time `split` is called. When comparing different estimators, this can
lead to overestimating the variance of the difference in performance between
the estimators::

    >>> from sklearn.naive_bayes import GaussianNB
    >>> from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.model_selection import KFold
    >>> from sklearn.model_selection import cross_val_score
    >>> import numpy as np

    >>> rng = np.random.RandomState(0)
    >>> X, y = make_classification(random_state=rng)
    >>> cv = KFold(shuffle=True, random_state=rng)
    >>> lda = LinearDiscriminantAnalysis()
    >>> nb = GaussianNB()

    >>> for est in (lda, nb):
    ...     print(cross_val_score(est, X, y, cv=cv))
    [0.8  0.75 0.75 0.7  0.85]
    [0.85 0.95 0.95 0.85 0.95]


Directly comparing the performance of the
:class:`~sklearn.discriminant_analysis.LinearDiscriminantAnalysis` estimator
vs the :class:`~sklearn.naive_bayes.GaussianNB` estimator on each fold would
be a mistake: **the splits on which the estimators are evaluated are
different**. Indeed, :func:`~sklearn.model_selection.cross_val_score` will
internally call `cv.split` on the same
:class:`~sklearn.model_selection.KFold` instance, but the splits will be
different each time. This is also true for any tool that performs model
selection via cross-validation, e.g.
:class:`~sklearn.model_selection.GridSearchCV` and
:class:`~sklearn.model_selection.RandomizedSearchCV`: scores are not
comparable fold-to-fold across different calls to `search.fit`, since
`cv.split` would have been called multiple times. Within a single call to
`search.fit`, however, fold-to-fold comparison is possible since the search
estimator only calls `cv.split` once.

For comparable fold-to-fold results in all scenarios, one should pass an
integer to the CV splitter: `cv = KFold(shuffle=True, random_state=0)`.

.. note::
    While fold-to-fold comparison is not advisable with `RandomState`
    instances, one can however expect that average scores allow to conclude
    whether one estimator is better than another, as long as enough folds and
    data are used.

.. note::
    What matters in this example is what was passed to
    :class:`~sklearn.model_selection.KFold`. Whether we pass a `RandomState`
    instance or an integer to :func:`~sklearn.datasets.make_classification`
    is not relevant for our illustration purpose. Also, neither
    :class:`~sklearn.discriminant_analysis.LinearDiscriminantAnalysis` nor
    :class:`~sklearn.naive_bayes.GaussianNB` are randomized estimators.

General recommendations
=======================

Getting reproducible results across multiple executions
.......................................................

In order to obtain reproducible (i.e. constant) results across multiple
*program executions*, we need to remove all uses of `random_state=None`, which
is the default. The recommended way is to declare a `rng` variable at the top
of the program, and pass it down to any object that accepts a `random_state`
parameter::

    >>> from sklearn.ensemble import RandomForestClassifier
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.model_selection import train_test_split
    >>> import numpy as np

    >>> rng = np.random.RandomState(0)
    >>> X, y = make_classification(random_state=rng)
    >>> rf = RandomForestClassifier(random_state=rng)
    >>> X_train, X_test, y_train, y_test = train_test_split(X, y,
    ...                                                     random_state=rng)
    >>> rf.fit(X_train, y_train).score(X_test, y_test)
    0.84

We are now guaranteed that the result of this script will always be 0.84, no
matter how many times we run it. Changing the global `rng` variable to a
different value should affect the results, as expected.

It is also possible to declare the `rng` variable as an integer. This may
however lead to less significant cross-validation results, as we will see in
the next section.

.. note::
    We do not recommend setting the global `numpy` seed by calling
    `np.random.seed(0)`. See `here
    <https://stackoverflow.com/questions/5836335/consistently-create-same-random-numpy-array/5837352#comment6712034_5837352>`_
    for a discussion.

Significance of cross-validation results
........................................

When we evaluate a randomized estimator performance by cross-validation, we
want to make sure that the estimator can yield accurate predictions for new
data, but we also want to make sure that the estimator is robust w.r.t. its
random initialization. For example, we would like the random weights
initialization of a :class:`~sklearn.linear_model.SGDCLassifier` to be
consistently good across all folds: otherwise, when we train that estimator
on new data, we might get unlucky and the random initialization may lead to
bad performance. Similarly, we want a random forest to be robust w.r.t the
set of randomly selected features that each tree will be using.

For these reasons, it is preferable to evaluate the cross-validation
preformance by letting the estimator use a different RNG on each fold. This
is done by passing a `RandomState` instance (or None) to the estimator
initialization.

When we pass an integer, the estimator will use the same RNG on each fold: if
we get a good (or bad) CV performance, it might just be because we got lucky
(or unlucky) with that specific seed. Passing instances leads to more
significant CV results, and makes the comparison between various algorithms
fairer. It also helps limiting the temptation to treat the estimator's RNG as
a hyper-parameter that can be tuned.

Whether we pass `RandomState` instances or integers to CV splitters has no
impact on significance, as long as `split` is only called once. When `split`
is called multiple times, fold-to-fold comparison isn't possible anymore. As
a result, passing integer to CV splitters is usually safer and covers most
use-cases.
