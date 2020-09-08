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
entry, and describes good practice and common pitfalls w.r.t. to this
subtle parameter.

The `random_state` parameter determines whether multiple calls to :term:`fit`
(for estimators) or to :term:`split` (for cv spliters) will produce the same
results, according to these rules:

- If an int is passed, calling `fit()` or `split()` multiple times yields the
  same results.
- If `None` or a `RandomState` instance is passed: `fit()` and `split()` will
  yield different results each time they are called, producing an event chain
  of maximal entropy. `None` is the default value for all `random_state`
  parameters.

  .. note::
      Since passing `random_state=None` is equivalent to passing the global
      `RandomState` instance from `numpy`
      (`random_state=np.random.mtrand._rand`), we will not explicitly mention
      `None` here, but everything that applies to instances also applies to
      using `None`.

.. warning:: TLDR

    Unless you know what you are doing, we strongly recommend to use integers
    as the `random_state` parameter of estimators and cv splitters. Leaving
    the default (`None`) or using `RandomState` instances is sometimes
    useful, but can have surprising effects.

Getting reproducible results across multiple executions
-------------------------------------------------------

In order to obtain reproducible (i.e. constant) results across multiple
program executions, the easiest way is to set a `rng` variable to an
integer value at the top of your program, and pass it down to any object that
accepts a `random_state` parameter::

    >>> from sklearn.ensemble import RandomForestClassifier
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.model_selection import train_test_split

    >>> rng = 0  # or any other integer
    >>> X, y = make_classification(random_state=rng)
    >>> rf = RandomForestClassifier(random_state=rng)
    >>> X_train, X_test, y_train, y_test = train_test_split(
    ...     X, y, random_state=rng)
    >>> rf.fit(X_train, y_train).score(X_test, y_test)
    0.92

We are now guaranteed that the result of this script will always be 0.92, no
matter how many times we run it. However, changing the global `rng` variable
should affect the results, as expected.

We could have used a `RandomState` instance instead of an integer::

    >>> import numpy as np
    >>> rng = np.random.RandomState(0)

This too would give us reproducible results across executions. However, using
`RandomState` instances may have surprising and unwanted effects that will be
described later.

Using None or `RandomState` instances
-------------------------------------

Estimators
..........

As described above, passing instances means that calling `fit()` multiple
times will not yield the same results, even if the estimator is fitted on the
same data, with the same hyper-parameters::

    >>> from sklearn.ensemble import RandomForestClassifier
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.model_selection import train_test_split
    >>> import numpy as np

    >>> rng = np.random.RandomState(0)
    >>> X, y = make_classification(random_state=rng)
    >>> rf = RandomForestClassifier(max_features=2, max_samples=10,
    ...                             random_state=rng)
    >>> X_train, X_test, y_train, y_test = train_test_split(
    ...     X, y, random_state=rng)
    >>> rf.fit(X_train, y_train).score(X_test, y_test)
    0.76
    >>> rf.fit(X_train, y_train).score(X_test, y_test)
    0.8

.. note::
    `score()` is not a random procedure. Only `fit()` has randomness.

We can see from the snippet above that `rf.fit()` has produced different
models, even if the data was the same. This is because the RNG of the
estimator is consumed when `fit()` is called, and this consumed (mutaded) RNG
will be used in the subsequent calls to `fit()`.

If we had passed an int to the `random_state` parameter of the
:class:`~sklearn.ensemble.RandomForestClassifier`, we would have obtained the
same models, and thus the same scores each time. When we pass an int, the
same RNG is used across all calls to `fit()`. What internally happens is that
even though the RNG is consumed when `fit()` is called, it is always reset to
its original state at the beginning of `fit()`.

CV splitters
............

Randomized cv splitters have a similar behavior when a `RandomState`
instance is passed::

    >>> from sklearn.model_selection import KFold
    >>> import numpy as np
    >>> X = np.arange(10)
    >>> rng = np.random.RandomState(0)
    >>> cv = KFold(n_splits=2, shuffle=True, random_state=rng)
    >>> for train, test in cv.split(X):
    ...     print(train, test)
    [0 3 5 6 7] [1 2 4 8 9]
    [1 2 4 8 9] [0 3 5 6 7]
    >>> for train, test in cv.split(X):
    ...     print(train, test)
    [0 4 6 7 8] [1 2 3 5 9]
    [1 2 3 5 9] [0 4 6 7 8]

We can see that the splits are different from the second time `split()` is
called. This may lead to wrong results if you compare the performance of
multiple estimators by calling `split()` many times: the estimators will not
be evaluated on the same folds, and performances will not be comparable.
Using an int is usually much safer.


Common pitfalls and subtleties
------------------------------

While the rules that govern the `random_state` parameter are seemingly simple,
they do however have some subtle implications. In some cases, this can even
lead to wrong conclusions.

Estimators
..........

**Differences in cross-validation procedures**

Depending on what is passed as the `random_state` parameter, estimators may
behave very differently, especially in cross-validation procedures. Consider
the following snippet::

    >>> from sklearn.ensemble import RandomForestClassifier
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.model_selection import cross_val_score
    >>> import numpy as np
    >>> 
    >>> X, y = make_classification(random_state=0)
    >>> rf_inst = RandomForestClassifier(random_state=np.random.RandomState(0))
    >>> cross_val_score(rf_inst, X, y)
    array([0.9 , 0.95, 0.95, 0.9 , 0.9 ])
    >>> rf_123 = RandomForestClassifier(random_state=123)
    >>> cross_val_score(rf_123, X, y)
    array([0.85, 0.95, 0.95, 0.9 , 0.9 ])

We see that the cross-validated scores of `rf_inst` and `rf_123` are
different, as should be expected since we didn't pass the same `random_state`
parameter. However, the difference between these scores is more subtle that
it may look, and **the cross-validation procedures that were performed by**
:func:`~sklearn.model_selection.cross_val_score` **significantly differ in
each case**:

- Since `rf_123` was passed an int, every call to `fit()` uses the same RNG:
  the same (random) subset of features will be used across all folds to fit
  the random forest.
- Since `rf_inst` was passed a `RandomState` instance, each call to `fit()`
  starts from a different RNG, and the randomly sampled subset of feature
  will be different for each of the 5 folds of the CV procedure.

Here, neither procedure is inherently wrong, and one might prefer one over
the other depending on the task at hand. It is however important to
understand how these procedures differ.

.. note::
    Here, :func:`~sklearn.model_selection.cross_val_score` will use a
    non-randomized cv splitter (as is the default), so both estimators will
    be evaluated on the same splits. Also, whether we pass an int or an
    instance to :func:`~sklearn.datasets.make_classification` isn't relevant
    for our illustration purpose: what matters is what we pass to the
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
will still be different models, even after calling `fit(X, y)` on the same
data. Moreover, `a` and `b` will influence each-other since they share the
same internal RNG: calling `a.fit()` will consume `b`'s RNG, and calling
`b.fit()` will consume `a`'s RNG, since they are the same.

If an int were passed, `a` and `b` would be exact clones and they would not
influence each other.

This is an important thing to remember because :func:`~sklearn.clone` is
called everywhere in scikit-learn tools: in particular, most meta-estimators
that accept non-fitted estimators will in fact call :func:`~sklearn.clone`
internally (:class:`~sklearn.model_selection.GridSearchCV`,
:class:`~sklearn.ensemble.StackingClassifier`,
:class:`~sklearn.calibration.CalibratedClassifierCV`, etc.).

CV splitters
............

When passed a `RandomState` instance, cv splitters yield different splits
each time `split()` is called. This can lead to dramatic mistakes when
comparing the performance of different estimators::

    >>> from sklearn.ensemble import RandomForestClassifier
    >>> from sklearn.ensemble import GradientBoostingClassifier
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.model_selection import KFold
    >>> from sklearn.model_selection import cross_val_score
    >>> import numpy as np
    >>> 
    >>> rng = np.random.RandomState(0)
    >>> X, y = make_classification(random_state=rng)
    >>> rf = RandomForestClassifier(random_state=rng)
    >>> gbdt = GradientBoostingClassifier(random_state=rng)
    >>> cv = KFold(shuffle=True, random_state=rng)
    >>> for est in (rf, gbdt):
    ...     print(cross_val_score(est, X, y, cv=cv))
    [0.85 0.95 0.9  0.95 0.95]
    [0.85 0.7  0.95 0.8  0.85]

Directly comparing the performance of the random forest vs the gradient
boosting estimator would be a methodological mistake: **the splits on which
the estimators are evaluated are different**. Indeed,
:func:`~sklearn.model_selection.cross_val_score` will internally call
`cv.split()` on the same :class:`~sklearn.model_selection.KFold` instance,
but the splits will be different each time. This is relevant for any tool
that performs model selection via cross-validation, including
:class:`~sklearn.model_selection.GridSearchCV` and
:class:`~sklearn.model_selection.RandomizedSearchCV`.

For comparable results, one should pass an int to
:class:`~sklearn.model_selection.KFold`: `KFold(shuffle=True,
random_state=0)`.

.. note::
    What matters in this example is what was passed to
    :class:`~sklearn.model_selection.KFold`. Whether we pass a `RandomState`
    instance or an int to :func:`~sklearn.datasets.make_classification` or to
    the estimators is not relevant for our illustration purpose. It does
    however have an impact on the cross-validation procedure as explained
    above, but it isn't what makes the comparison incorect.
