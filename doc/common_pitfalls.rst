Common Pitfalls
===============

The purpose of this guide is to illustrate some common pitfalls and
anti-patterns that occur when using ``scikit-learn``. It provides
examples of what **not** to do, along with a corresponding correct
example.

Inconsistent preprocessing
--------------------------

scikit-learn provides a library of :ref:`data-transforms`, which
may clean (see :ref:`preprocessing`), reduce
(see :ref:`data_reduction`), expand (see :ref:`kernel_approximation`)
or generate (see :ref:`feature_extraction`) feature representations.
If these data transforms are used when training a model, they also
must be used on subsequent datasets, whether it's test data or
data in a production system. Otherwise, the feature space will change,
and the model will not be able to perform effectively.

For the following example, let's create a synthetic dataset with a
single feature::

    >>> from sklearn.metrics import mean_squared_error
    >>> from sklearn.datasets import make_regression
    >>> from sklearn.model_selection import train_test_split
    >>> from sklearn.pipeline import make_pipeline
    >>> from sklearn.preprocessing import StandardScaler
    >>> from sklearn.linear_model import LinearRegression
    ...
    >>> random_state = 42
    >>> X, y = make_regression(random_state=random_state, n_features=1, noise=1)
    >>> X_train, X_test, y_train, y_test = train_test_split(
    ... X, y, test_size=0.4, random_state=random_state)

**Wrong**

The train dataset is scaled, but not the test dataset, so model
performance on the test dataset is worse than expected::

    >>> scaler = StandardScaler()
    >>> X_train_transformed = scaler.fit_transform(X_train)
    >>> model = LinearRegression().fit(X_train_transformed, y_train)
    >>> mean_squared_error(y_test, model.predict(X_test))
    62.80...

**Right**

A :class:`Pipeline <sklearn.pipeline.Pipeline>` makes it easier to chain
transformations with estimators, and decreases the possibility of
forgetting a transformation::

    >>> model = make_pipeline(StandardScaler(), LinearRegression()).fit(
    ... X_train, y_train)
    >>> mean_squared_error(y_test, model.predict(X_test))
    0.90...

Next steps
----------

We have briefly covered some common pitfalls and anti-patterns, but
there is much more to ``scikit-learn``!

Please refer to our :ref:`user_guide` for details on all the tools that we
provide. You can also find an exhaustive list of the public API in the
:ref:`api_ref`.

You can also look at our numerous :ref:`examples <general_examples>` that
illustrate the use of ``scikit-learn`` in many different contexts.

The :ref:`tutorials <tutorial_menu>` also contain additional learning
resources.
