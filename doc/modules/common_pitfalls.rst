==============================
Common methodological pitfalls
==============================


Feature Selection
=================

**Wrong Practice: Selecting features using Whole Dataset**

    >>> import numpy
    >>> from sklearn.feature_selection import SelectKBest
    >>> from sklearn.model_selection import cross_val_score
    >>> from sklearn.ensemble import GradientBoostingClassifier
    >>> n_samples, n_features, n_classes = 200, 10000, 2
    >>> X = numpy.random.randn(n_samples, n_features)
    >>> y = numpy.random.choice(n_classes, n_samples)
    >>> X_selected = SelectKBest(k=25).fit_transform(X, y)
    >>> print("Mean Accuracy without Pipeline:{0.2f}".format(cross_val_score(GradientBoostingClassifier
    ...      (), X_selected, y, cv=5).mean()*100))
    Mean Accuracy without Pipeline:73.5
    >>> cv_scores = cross_val_score(GradientBoostingClassifier(), X_selected, y, cv=5)
    >>> print("Cross-validation accuracy: {0.2f}+/-{0.2f}".format(cv_scores.mean(), cv_scores.std())
    Cross-validation accuracy: 0.73 +/- 0.04

Since the y class values are independently sampled from the input features
X, one should be surprised to see the GradientBoostingClassifier model
performing significantly better than the chance level.

The cross-validation procedure applied here is not applied to a model
trained on the original data X but instead reflects the performance of the GradientBoostingClassifier
trained on X_selected which is now dependent on y because of the supervised feature selection step. 
This explains why we get a good performance accuracy instead of chance level.

However what we are really interested in is quantifying the performance
of the joint model that combines both steps: the feature selection step
and the classification step. To do so we need to wrap the feature
selection step inside the cross-validation loop.

**Right Practice: Selecting features using Training Dataset**

    >>> from sklearn.pipeline import make_pipeline
    >>> pp = make_pipeline(SelectKBest(k=25), GradientBoostingClassifier())
    >>> print("Mean Accuracy with Pipeline:{0.2f}".format(cross_val_score(pp,
    ...       X, y, cv=5).mean()*100))
    Mean Accuracy with Pipeline:58.5
    >>> cv_scores = cross_val_score(pp, X, y, cv=5)
    >>> print("Cross-validation accuracy: {0.2f} +/- {0.2f}".format(
    ...     cv_scores.mean(), cv_scores.std()))
    Cross-validation accuracy: 0.57 +/- 0.04

Inside the pipeline, the features are selected using training set of each
cross-validation split. This results in more reliable accuracy of the
model. This accuracy is a better representative of the joint model that
includes both the feature selection step and the classifier.