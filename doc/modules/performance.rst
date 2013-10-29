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
linear algebra libraries optimizations etc.). Here we see on a setting
with few features that independently of estimator choice the bulk mode is
always faster by 2 orders of magnitude:

.. |atomic_prediction_speed| image::  ../auto_examples/applications/images/plot_prediction_latency_1.png
    :target: ../auto_examples/applications/plot_prediction_latency.html
    :scale: 80

.. centered:: |atomic_prediction_speed|

.. |bulk_prediction_speed| image::  ../auto_examples/applications/images/plot_prediction_latency_2.png
    :target: ../auto_examples/applications/plot_prediction_latency.html
    :scale: 80

.. centered:: |bulk_prediction_speed|

To benchmark different estimators for your case you can simply change the
`n_feats` parameter according to your case in this example:
:ref:`example_applications_plot_prediction_latency.py`. This should give you
an estimate of the order of magnitude of the prediction speed for your case.

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

In many real world applications the feature extraction process (i.e. turning
raw data like database rows or network packets into numpy arrays) governs the
overall prediction time. For example here on the Reuters text classification
task the vectorization that includes parsing SGML files, tokenizing the text
and hashing it into a common vector space is taking 5 to 30 times more time
than the actual prediction code, depending on the chosen model.

 .. |computation_time| image::  ../auto_examples/applications/images/plot_out_of_core_classification_3.png
    :target: ../auto_examples/applications/plot_out_of_core_classification.html
    :scale: 80

.. centered:: |computation_time|

In many cases it is thus recommended to carefully time and profile your
feature extraction code as it may be a good place to start optimizing when
your overall speed is too slow for your application. If needed,
you can consider rewriting the feature extraction part in a lower-level,
compiled language to further speed up the overall process. The fact that
most scikit-learn models are implemented using Cython and optimized,
compiled computing libraries under the hood make them usually pretty fast.
So optimizing the feature extraction step while keeping the prediction in
python with scikit-learn estimators is usually a good way to go as it allows
for easy experimentation on the modeling side without sacrificing performance.

Tips and Tricks
===============

Linear algebra libraries
------------------------
Compile Atlas for your architecture etc.

Model Compression
-----------------
`sparsify()` trick etc

