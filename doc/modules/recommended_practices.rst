=====================
Recommended Practices
=====================


Feature Selection
=================

**Wrong Practice: Selecting features using Whole Dataset**

    >>> X = numpy.random.randn(500,10000)
    >>> y = numpy.random.choice(2,500)
    >>> X_ = SelectKBest(k=25).fit_transform(X,y)
    >>> print("Mean Accuracy without Pipeline:{0.2f}".format(cross_val_score
    ...       (GradientBoostingClassifier(),X_,y,cv=3,n_jobs=8).mean()*100))
    Mean Accuracy without Pipeline: 64.79

The accuracy of model is almost 65%, 15% more than what would be expected from
random guessing.

**Right Practice: Selecting features using Training Dataset**

    >>> from sklearn.pipeline import make_pipeline
    >>> pp = make_pipeline(SelectKBest(k=25), GradientBoostingClassifier())
    >>> print("Mean Accuracy with Pipeline:{0.2f}".format(cross_val_score(pp,
    ...       X, y, cv=5).mean()*100))
    Mean Accuracy with Pipeline:44.77

Inside the pipeline, the features are selected using Training dataset. This
results in more realistic and reliable Accuracy of the model.
