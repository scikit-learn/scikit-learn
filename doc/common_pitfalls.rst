.. Places parent toc into the sidebar

:parenttoc: True

.. include:: includes/big_toc_css.rst

.. _common_pitfalls:

===============
Common pitfalls
===============

The purpose of this chapter is to illustrate some common pitfalls and
anti-patterns that occur when using scikit-learn. It provides
examples of what **not** to do, along with a corresponding correct
example.

Inconsistent preprocessing
==========================

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


    >>> from sklearn.datasets import make_regression
    >>> from sklearn.model_selection import train_test_split
    ...
    >>> random_state = 42
    >>> X, y = make_regression(random_state=random_state, n_features=1, noise=1)
    >>> X_train, X_test, y_train, y_test = train_test_split(
    ... X, y, test_size=0.4, random_state=random_state)

**Wrong**

The train dataset is scaled, but not the test dataset, so model
performance on the test dataset is worse than expected::

    >>> from sklearn.metrics import mean_squared_error
    >>> from sklearn.linear_model import LinearRegression
    >>> from sklearn.preprocessing import StandardScaler
    ...
    >>> scaler = StandardScaler()
    >>> X_train_transformed = scaler.fit_transform(X_train)
    >>> model = LinearRegression().fit(X_train_transformed, y_train)
    >>> mean_squared_error(y_test, model.predict(X_test))
    62.80...

**Right**

A :class:`Pipeline <sklearn.pipeline.Pipeline>` makes it easier to chain
transformations with estimators, and reduces the possibility of
forgetting a transformation::

    >>> from sklearn.pipeline import make_pipeline
    ...
    >>> model = make_pipeline(StandardScaler(), LinearRegression())
    >>> model.fit(X_train, y_train)
    Pipeline(steps=[('standardscaler', StandardScaler()),
                    ('linearregression', LinearRegression())])
    >>> mean_squared_error(y_test, model.predict(X_test))
    0.90...

.. _data_leakage:

Data leakage
============

Data leakage occurs when information that would not be available at prediction
time is used when building the model. This results in overly optimistic
performance estimates, for example from :ref:`cross-validation
<cross_validation>`, and thus poorer performance when the model is used
on actually novel data, for example during production.

A common cause is not keeping the test and train data subsets separate.
Test data should never be used to make choices about the model.
The general rule is to never call `fit` on test data. This is particularly
important for 'supervised transformations' where the `fit` method requires
both the data, `X` and the target values, `y`. While this may
sound obvious, this is easy to miss in some cases, for example when applying
certain pre-processing steps.

Although both train and test data subsets should receive the same preprocessing
transformation, it is important that these transformations are only learnt
from the training data. For example, if you have a
normalization step where you divide by the average value, the average should
be the average of the train subset, **not** the average of all the data. If the
test subset is included in the average calculation, information from the test
subset is influencing the model.

Including the test data when :ref:`tuning model hyperparameters <grid_search>`
will also inadvertently introduce information from the test data into the
model. Practically, this means that only the train data subset should be fed
into the `fit` method of :ref:`hyper_parameter_optimizers`, as is the case
with all scikit-learn estimators.

An example of data leakage during preprocessing is detailed below.

.. topic:: See Also:

  * :ref:`sphx_glr_auto_examples_model_selection_plot_nested_cross_validation_iris.py`

Data leakage during feature selection
-------------------------------------

A number of :ref:`feature_selection` functions are available in scikit-learn.
They can help remove irrelevant, redundant and noisy features as well as
improve your model build time and performance. As with any other type of
preprocessing, feature selection should **only** use the training data.
Including the test data in feature selection will optimistically bias your
model.

To demonstrate we will create this binary classification problem with
10,000 randomly generated features::

    >>> import numpy as np
    >>> n_samples, n_features, n_classes = 200, 10000, 2
    >>> rng = np.random.RandomState(42)
    >>> X = rng.standard_normal((n_samples, n_features))
    >>> y = rng.choice(n_classes, n_samples)

**Wrong**

Using all the data to perform feature selection results in an accuracy score
much higher than chance, even though our targets are completely random.
This randomness means that our `X` and `y` are independent and we thus expect
the accuracy to be around 0.5. However, since the feature selection step
'sees' the test data, the model has an unfair advantage. In the incorrect
example below we first use all the data for feature selection and then split
the data into training and test subsets for model fitting. The result is a
much higher than expected accuracy score::

    >>> from sklearn.model_selection import train_test_split
    >>> from sklearn.feature_selection import SelectKBest
    >>> from sklearn.ensemble import GradientBoostingClassifier
    >>> from sklearn.metrics import accuracy_score
    >>> X_selected = SelectKBest(k=25).fit_transform(X, y)
    >>> X_train, X_test, y_train, y_test = train_test_split(
    ...     X_selected, y, random_state=42)
    >>> gbc = GradientBoostingClassifier(random_state=1)
    >>> gbc.fit(X_train, y_train)
    GradientBoostingClassifier(random_state=1)
    >>> y_pred = gbc.predict(X_test)
    >>> score = accuracy_score(y_test, y_pred)
    >>> print(f"Accuracy: {score}")
    Accuracy: 0.76

**Right**

To prevent data leakage, it is good practice to split your data into train
and test subsets **first**. Feature selection can then be formed using just
the train dataset. Notice that whenever we use `fit` or `fit_transform`, we
only use the train dataset. The score is now what we would expect for the
data, close to chance::

    >>> X_train, X_test, y_train, y_test = train_test_split(
    ...     X, y, random_state=42)
    >>> select = SelectKBest(k=25)
    >>> X_train_selected = select.fit_transform(X_train, y_train)
    >>> gbc = GradientBoostingClassifier(random_state=1)
    >>> gbc.fit(X_train_selected, y_train)
    GradientBoostingClassifier(random_state=1)
    >>> X_test_selected = select.transform(X_test)
    >>> y_pred = gbc.predict(X_test_selected)
    >>> score = accuracy_score(y_test, y_pred)
    >>> print(f"Accuracy: {score}")
    Accuracy: 0.46

Another way to prevent data leakage is to use the
:class:`~sklearn.pipeline.Pipeline` to chain together the feature selection
and model estimators. The pipeline ensures that only the training data is
used when performing `fit` and the test data is used only for calculating the
accuracy score::

    >>> from sklearn.pipeline import make_pipeline
    >>> X_train, X_test, y_train, y_test = train_test_split(
    ...     X, y, random_state=42)
    >>> pipeline = make_pipeline(SelectKBest(k=25),
    ...                          GradientBoostingClassifier(random_state=1))
    >>> pipeline.fit(X_train, y_train)
    Pipeline(steps=[('selectkbest', SelectKBest(k=25)),
                    ('gradientboostingclassifier',
                    GradientBoostingClassifier(random_state=1))])
    >>> y_pred = pipeline.predict(X_test)
    >>> score = accuracy_score(y_test, y_pred)
    >>> print(f"Accuracy: {score.mean():.2f}")
    Accuracy: 0.46

The pipeline can also be fed into a cross-validation
function such as :func:`~sklearn.model_selection.cross_val_score`.
Again, the pipeline ensures that the correct data subset and estimator
method is used during fitting and predicting::

    >>> from sklearn.model_selection import cross_val_score
    >>> scores = cross_val_score(pipeline, X, y)
    >>> print(f"Mean accuracy: {scores.mean():.2f}+/-{scores.std():.2f}")
    Mean accuracy: 0.45+/-0.07

How to avoid data leakage
-------------------------

Below are some tips on avoiding data leakage:

* Always split the data into train and test subsets first, particularly
  before any preprocessing steps.
* Never include test data when using estimator `fit` and `fit_transform`
  methods. Using all the data, e.g., `fit(X)`, can result in overly optimistic
  scores. Only using the test data, e.g., `fit(X_test)`, can harm performance
  as preprocessing or model fitting is only performed using the, generally
  smaller, test subset.
  Conversely, the `transform` method should be used on both train and test
  subsets as the same preprocessing should be applied to all the data.
  This can be achieved by using `fit_transform`, which combines the `fit` and
  `transform` methods, on the train subset and `transform` on the test
  subset.
* The scikit-learn :ref:`pipeline <pipeline>` is a great way to prevent data
  leakage as it ensures that the appropriate method is performed on the
  correct data subset. The pipeline is ideal for use in cross-validation
  and hyper-parameter tuning functions.
