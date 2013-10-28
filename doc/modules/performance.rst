.. _performance:

===========
Performance
===========

For some applications the performance (mainly speed and throughput) of
estimators is crucial. We will review here the orders of magnitude you can
expect from a number of scikit-learn estimators in different contexts and
provide some tips and tricks for overcoming performance bottlenecks.

Prediction Speed
================

One of the most straight-forward concerns one may have when using/choosing a
machine learning toolkit is the speed at which predictions can be made in a
production environment.

Scikit-learn provides a number of benchmark utilities :ref:`benchmarks link
 <bench_link>` and :ref:`example link <ex_link>` that can be used to estimate
the order of magnitude you can expect.

The main factors that influence the prediction speed are:
  1. Number of features
  2. Matrix type
  3. Model complexity
A last major parameter is also the possibility to do predictions in bulk or
one-at-a-time mode.

Bulk versus Atomic mode
-----------------------
In general doing predictions in bulk (many instances at the same time) is
more efficient for a number of reasons (branching predictability, CPU cache,
linear algebra libraries optimizations etc.). Here we see that independently
of estimator choice the bulk mode is always faster by 2 orders of magnitude:

.. |atomic_prediction_speed| image::  ../auto_examples/applications/images/plot_prediction_latency_1.png
    :target: ../auto_examples/applications/plot_prediction_latency.html
    :scale: 80

.. centered:: |atomic_prediction_speed|

.. |bulk_prediction_speed| image::  .
./auto_examples/applications/images/plot_prediction_latency_2.png
    :target: ../auto_examples/applications/plot_prediction_latency.html
    :scale: 80

.. centered:: |bulk_prediction_speed|

Influence of Dimensionality
---------------------------
tbd

Influence of the number of Features
-----------------------------------
tbd

Influence of Matrix Type
------------------------
tbd (CSR vs dense vs ...)

Prediction Throughput
=====================

Feature Extraction Speed
========================

Tips and Tricks
===============

Linear algebra libraries
------------------------

