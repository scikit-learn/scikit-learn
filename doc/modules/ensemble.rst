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

    **Examples:** Bagging methods, :ref:`forests of randomized trees <forest>`, ...

  - By contrast, in *boosting methods*, models are built sequentially and one
    tries to reduce the bias of the combined model. The motivation is to combine
    several weak models to prodce a powerful ensemble.

    **Examples:** AdaBoost, Least Squares Boosting, Gradient Tree Boosting, ...


.. _forest:

Forests of randomized trees
===========================

TODO
