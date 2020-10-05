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

.. _data_leakage:

Data leakage
============

Data leakage occurs when information that would not be available at prediction
time is used when building the model. This results in overly optimsitic
performance estimates, for example from :ref:`cross-validation
<cross_validation>`, and thus poorer performance when the model is used
on actually novel data, for example during production.

A common cause is not keeping the test and train data subsets separate.
Test data should never be used to make choices about the model.
The general rule is to never call `fit` on testing data. This is particularly
important for 'supervised transformations' where the `fit` method requires
both the data, `X` and the target values, `y`. While this may
sound obvious, this is easy to miss in some cases, for example when applying
certain pre-processing steps.

Although both train and test data subsets should receive the same preprocessing
transformation, it is important that these transformations are only learnt
from the training data. For example, if you have a
normalization step where you divide by the average value, the average should
be the average of the train subset, **not** the average of all the data. If the
test subset was included in the average calculation, information from the test
subset is influencing the model.

Including the test data when finding the best model hyperparameters will
also inadvertantly introduce information from the test data into the model.
Practically, this means that only the train data subset should be fed into
the `fit` method of :ref:`hyper_parameter_optimizers`, as is the case with all
scikit-learn estimators.

Some examples of common data leakage pitfalls are detailed below.

Data leakage during feature selection
-------------------------------------

A number of :ref:`feature_selection` functions are available in scikit-learn.
They can help remove irrelevant, redundant and noisy features as well as
improve your model build time and performance. As with any other type of
preprocessing, feature selection should
**only** use the training data. Including the test data in feature selection
will optimistically bias your model.

To demonstrate we will create this binary classification problem with
1000 randomly generated features::

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
'sees' the test data, the model has an unfair advantage. Below we use all the
data for feature selection before splitting the data into training and test
subsets for model fitting. The result is a much higher than expected
accuracy score::

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
and model estimators. This pipeline can then be fed into a cross-validation
function such as :func:`~sklearn.model_selection.cross_val_score`. The
pipeline ensures that only the training data when performing `fit` and
the test data will only be used for calculating the accuracy score.
Using :func:`~sklearn.model_selection.cross_val_score` also provides us
information about the variance of the accuracy score::

    >>> from sklearn.pipeline import make_pipeline
    >>> from sklearn.model_selection import cross_val_score
    >>> pipeline = make_pipeline(SelectKBest(k=25),
    ...                          GradientBoostingClassifier(random_state=1))
    >>> scores = cross_val_score(pipeline, X, y)
    >>> print(f"Mean Accuracy: {scores.mean():.2f}+/-{scores.std():.2f}")
    Mean Accuracy: 0.45+/-0.07

Data leakage during imputation
------------------------------

There are a number of methods to impute missing values in data. For example,
:class:`~sklearn.impute.KNNIMputer` uses the mean value from neighbors to
impute missing values. Only the train data should be used to calculate this
mean value, as including the test data in the mean calculation will introduce
information about the test data into the model.

To demonstrate this, we will use the :ref:`diabetes_dataset` and
artificially introduce missing values. The smallest 100 `y` values are 10 times
more likely to be missing, simulating 'missing not at random'::

    >>> import numpy as np
    >>> from sklearn.datasets import load_diabetes
    >>> X, y = load_diabetes(return_X_y=True)
    >>> n_samples, n_features = X.shape
    >>> indx = np.argsort(y)
    >>> X_sorted = X[indx, :]
    >>> y_sorted = y[indx]
    >>> rng = np.random.RandomState(42)
    >>> mask1 = rng.binomial(n=1, p=0.1, size=(100, n_features))
    >>> mask2 = rng.binomial(n=1, p=0.01, size=(n_samples-100, n_features))
    >>> full_mask = np.vstack((mask1, mask2)).astype(bool)
    >>> X_missing = X_sorted.copy()
    >>> X_missing[full_mask] = np.nan

**Wrong**

Using all the data to calculate impute the missing values, results in an
overly optimsitic :math:`R^2`::

    >>> from sklearn.impute import KNNImputer
    >>> from sklearn.model_selection import train_test_split
    >>> from sklearn.ensemble import GradientBoostingRegressor
    >>> from sklearn.metrics import r2_score
    >>> X_impute = KNNImputer().fit_transform(X_missing)
    >>> X_train, X_test, y_train, y_test = train_test_split(
    ...     X_impute, y, random_state=42)
    >>> gbr = GradientBoostingRegressor(random_state=1)
    >>> gbr.fit(X_train, y_train)
    GradientBoostingRegressor(random_state=1)
    >>> y_pred = gbr.predict(X_test)
    >>> score = r2_score(y_test, y_pred)
    >>> print(f"R2 score: {score:.3f}")
    R2 score: -0.188

**Right**

As above, splitting your data into test and train subsets should be done
first. This enables imputation to be performed using just the train subset,
which can then used to fit our model. The :math:`R^2` score is now much
less accurate::

    >>> X_train, X_test, y_train, y_test = train_test_split(
    ...     X_missing, y, random_state=42)
    >>> impute = KNNImputer()
    >>> X_train_impute = impute.fit_transform(X_train)
    >>> gbr = GradientBoostingRegressor(random_state=1)
    >>> gbr.fit(X_train_impute, y_train)
    GradientBoostingRegressor(random_state=1)
    >>> X_test_impute = impute.transform(X_test)
    >>> y_pred = gbr.predict(X_test_impute)
    >>> score = r2_score(y_test, y_pred)
    >>> print(f"R2 score: {score:.3f}")
    R2 score: -0.159

The :class:`~sklearn.pipeline.Pipeline` is another way to prevent data
leakage. It chains together the imputation and model estimators and ensures
that the correct data subset is used for fit, transform and predict when
used for example, in a cross-validation function. This is shown below
along with the mean and standard deviation of the :math:`R^2` scores from
cross-validation::

    >>> from sklearn.pipeline import make_pipeline
    >>> from sklearn.model_selection import cross_val_score
    >>> pipeline = make_pipeline(KNNImputer(),
    ...                          GradientBoostingRegressor(random_state=1))
    >>> scores = cross_val_score(pipeline, X_missing, y)
    >>> print(f"Mean R2: {scores.mean():.3f}+/-{scores.std():.2f}")
    Mean R2: -0.220+/-0.09

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
  subsets as the same proprocessing should be applied to all the data.
  This can be achieved by using `fit_transform`, which combines the `fit` and
  `transform` methods, on the train subset and `transform` on the test
  subset.
* The scikit-learn :ref:`pipeline <pipeline>` is a great way to prevent data
  leakage as it ensures that the appropriate method is performed on the
  correct data subset. The pipeline is ideal for use in cross-validation
  and hyper-parameter tuning functions.
