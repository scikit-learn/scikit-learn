.. _wrong_feature_selection:

================================
Feature selection using all data
================================

A number of :ref:`feature_selection` functions are available in scikit-learn.
They can help remove irrelevant, redundant and noisy features as well as
improve your model build time and performance. Feature selection should
**only** use the training data. Indeed, test data should never be used to make
choices about the model. Including the test data in feature selection will
optimistically bias your model.

To demonstrate we will create a binary classification problem with
1000 randomly generated features::

    >>> from numpy.random import default_rng
    >>> n_samples, n_features, n_classes = 200, 10000, 2
    >>> rng = default_rng(42)
    >>> X = rng.standard_normal((n_samples, n_features))
    >>> y = rng.choice(n_classes, n_samples)

**Wrong**

Using all the data to perform feature selection results in an accuracy score
much higher than chance. Since our `X` and `y` are independent, we expect
performance to be around 0.5. However, since the feature selection step
'sees' the test data, the model has an unfair advantange.

    >>> from sklearn.feature_selection import SelectKBest
    >>> from sklearn.model_selection import cross_val_score
    >>> from sklearn.ensemble import GradientBoostingClassifier
    >>> X_selected = SelectKBest(k=25).fit_transform(X, y)
    >>> scores = cross_val_score(GradientBoostingClassifier(random_state=1),
    ...                          X_selected, y, cv=5)
    >>> print(f"Mean cross-validation accuracy: {scores.mean()}")
    Mean cross-validation accuracy: 0.775

**Right**

To ensure that feature selection is performed using only the train dataset
we will use a :class:`~sklearn.pipeline.Pipeline` to chain together the
feature selection and model. Feeding our pipeline into
:func:`~sklearn.model_selection.cross_val_score` ensures that only the
training data is used for feature selection (and when fitting our model).
The test data will only be used for calculating the accuracy score. The
score is now what we would expect for the data, close to chance::

    >>> from sklearn.pipeline import make_pipeline
    >>> pipeline = make_pipeline(SelectKBest(k=25),
    ...                          GradientBoostingClassifier(random_state=1))
    >>> scores = cross_val_score(pipeline, X, y, cv=5)
    >>> print(f"Mean Accuracy: {scores.mean()}")
    Mean Accuracy: 0.495
