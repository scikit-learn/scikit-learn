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
    combined model is hence usually better than any of the single model because
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
specifically designed for trees.

In random forests (see :class:`RandomForestClassifier` and
:class:`RandomForestRegressor` classes), each tree in the ensemble is built
on a bootstrap sample drawn from the training set. In addition, when splitting
a node during the construction of the tree, the split that is chosen is no
longer the best split among all features. Instead, the split that is picked
is the best split among a random subset of the features.

In extra-trees...

.. topic:: References

 * Leo Breiman, "Random Forests", Machine Learning, 45(1), 5-32, 2001.

 * Pierre Geurts, Damien Ernst., and Louis Wehenkel, "Extremely randomized
   trees", Machine Learning, 63(1), 3-42, 2006.
