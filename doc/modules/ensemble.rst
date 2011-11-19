.. _ensemble:

================
Ensemble methods
================

.. currentmodule:: sklearn.ensemble

The goal of **ensemble methods** is to combine the predictions of several
models built with a given learning algorithm in order to improve with respect
to the use of a single model.

Two families of ensemble methods are usually distinguished:

  - In *averaging methods*, the driving principle is to build several models
    independently and then to average their predictions. On average, the
    combined model is usually better than any of the single model because
    its variance is reduced.

    **Examples:** Bagging methods, :ref:`Forests of randomized trees <forest>`, ...

  - By contrast, in *boosting methods*, models are built sequentially and one
    tries to reduce the bias of the combined model. The motivation is to combine
    several weak models to produce a powerful ensemble.

    **Examples:** AdaBoost, Least Squares Boosting, Gradient Tree Boosting, ...


.. _forest:

Forests of randomized trees
===========================

The ``sklearn.ensemble`` module includes two averaging algorithms based on
randomized :ref:`decision trees <tree>`: the RandomForest algorithm and the
Extra-Trees method. Both algorithms are perturb-and-combine techniques
specifically designed for trees::

    >>> from sklearn.ensemble import RandomForestClassifier
    >>> X = [[0, 0], [1, 1]]
    >>> Y = [0, 1]
    >>> clf = RandomForestClassifier(n_estimators=10)
    >>> clf = clf.fit(X, Y)

In random forests (see :class:`RandomForestClassifier` and
:class:`RandomForestRegressor` classes), each tree in the ensemble is built from
a bootstrap sample drawn from the training set. In addition, when splitting a
node during the construction of the tree, the split that is chosen is no longer
the best split among all features. Instead, the split that is picked is the best
split among a random subset of the features. As a result of this randomness, the
bias of the forest usually slightly increases (with respect to the bias of a
single non-random tree) but, due to averaging, its variance also decreases,
usually more than compensating for the increase in bias, hence yielding an
overall better model.

In extra-trees (see :class:`ExtraTreesClassifier` and
:class:`ExtraTreesRegressor` classes), randomness goes one step further in the
way splits are computed. As in random forests, a random subset of candidate
features is used, but instead of looking for the most discriminative thresholds,
thresholds are drawn at random for each candidate feature and the best of these
randomly-generated thresholds is picked as the splitting rule. This usually
allows to reduce the variance of the model a bit more, at the expense of a
slightly greater increase in bias::

    >>> from sklearn.cross_validation import cross_val_score
    >>> from sklearn.datasets import make_blobs
    >>> from sklearn.ensemble import RandomForestClassifier
    >>> from sklearn.ensemble import ExtraTreesClassifier
    >>> from sklearn.tree import DecisionTreeClassifier

    >>> X, y = make_blobs(n_samples=10000, n_features=10, centers=100,
    ...     random_state=0)

    >>> clf = DecisionTreeClassifier(max_depth=None, min_split=1,
    ...     random_state=0)
    >>> scores = cross_val_score(clf, X, y)
    >>> scores.mean()                             # doctest: +ELLIPSIS
    0.978...

    >>> clf = RandomForestClassifier(n_estimators=10, max_depth=None,
    ...     min_split=1, random_state=0)
    >>> scores = cross_val_score(clf, X, y)
    >>> scores.mean()                             # doctest: +ELLIPSIS
    0.992...

    >>> clf = ExtraTreesClassifier(n_estimators=10, max_depth=None,
    ...     min_split=1, random_state=0)
    >>> scores = cross_val_score(clf, X, y)
    >>> scores.mean()                             # doctest: +ELLIPSIS
    0.999...

The main parameters to adjust when using these methods is ``n_estimators``
and ``max_features``. The former is the number of trees in the
forest. The  larger the better, but also the longer it will take to
compute. The latter is the size of the random subsets of features to
consider when splitting a node. The lower the greater the reduction of
variance, but also the greater the increase in bias. Empiricial good
default values are ``max_features=n_features`` in random forests, and
``max_features=sqrt(n_features)`` in extra-trees (where ``n_features``
is the number of features in the data). The best results are also
usually reached when setting ``max_depth=None`` in combination with
``min_split=1`` (i.e., when fully developping the trees).

Finally, note that bootstrap samples are used by default in random forests
(``bootstrap=True``) while the default strategy is to use the original
datasets for building extra-trees (``bootstrap=False``).

.. topic:: Examples:

 * :ref:`example_ensemble_plot_forest_iris.py`

.. topic:: References

 * Leo Breiman, "Random Forests", Machine Learning, 45(1), 5-32, 2001.

 * Pierre Geurts, Damien Ernst., and Louis Wehenkel, "Extremely randomized
   trees", Machine Learning, 63(1), 3-42, 2006.
