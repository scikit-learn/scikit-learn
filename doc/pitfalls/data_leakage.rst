.. _data_leakage:

============
Data leakage
============

Data leakage occurs when information that would not be available at prediction
time, is used when building the model. This results in overly optimstic
performance estimates, for example from :ref:`cross_validation`, but
poorer performance when the model is asked to predict on actually novel data,
for example during production.

A common cause is not keeping the test and train data subsets separate. Test
data should never be used to make choices about the model but this may
accidently occur during :ref:`preprocessing` or :ref:`grid_search`.

Although both train and test data subsets should undergo the same preprocessing
transformation, it is important that stateful transformations use only the
training subset to determine the 'state'. For example, if you have a
normalization step whereby you divide by the average value, the average should
be the average of the train subset, **not** the average of all the data. If the
test subset was included in the average calculation, information from the test
subset is influencing the model. Other types of preprocessing such as
:ref:`impute` and :ref:`polynomial_features` should also only utilize train
data.

Including the test data when finding the best model hyperparameters will
also inadvertantly introduce information from the test data into the model.
Practically, this means that only the train data subset should be fed into the
`fit` method of ref:`hyper_parameter_optimizers`.

Examples of common data leakage pitfalls are detailed below.

Data leakage during feature selection
=====================================

A number of :ref:`feature_selection` functions are available in scikit-learn.
They can help remove irrelevant, redundant and noisy features as well as
improve your model build time and performance. Feature selection should
**only** use the training data. Including the test data in feature selection
will optimistically bias your model.

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


