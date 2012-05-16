.. _linear_gaussian:

=====================
Linear-Gaussian Model
=====================

.. currentmodule:: sklearn.linear_gaussian

`sklearn.linear_gaussian` is the continuous domain sister of the methods implemented in `sklearn.hmm`.  Like the Hidden Markov Model, The Linear-Gaussian Model is a generative probability model for explaining a sequence of measurements observed one after another.  The model assumes that there exists some unobserved "true state" of the system, and that these measurements are simply noisy version of that true state.

In particular, the Linear-Gaussian Model and Hidden Markov Model both assume that the true state of the system, :math:`x_t`, and the measurements, :math:`z_t`, are influenced by each other as described in the following diagram.

.. figure:: ../images/hmm.jpg
  :align center
  :width: 600 px

Where these two models differ is in what these variables represent.  While in a Hidden Markov Model both :math:`x_t` and :math:`z_t` are one of a finite set of values, in a Linear-Gaussian Model both are vectors of real numbers.  

The fundamental problems in a Linear-Gaussian Model are identical to those in the Hidden Markov Model, namely:

* Given the model parameters and a sequence of measurements, estimate the most likely sequence of hidden states
* Given the model parameters and a sequence of measurements, estimate the likelihood of the measurements
* Given the measurements, estimate the model parameters.

This submodule implements the Kalman Filter and Kalman Smoother, two algorithms for solving the first two goals.  The third goal is solved by the EM algorithm as applied to the Linear-Gaussian model.

Mathematical Formulation
========================

In order to understand when the algorithms in this module will be effective, it is important to understand what assumptions are being made.  In words, the Linear-Gaussian model assumes that for all :math:`t = 0, \ldots, T`,

* :math:`x_0` is distributed according to a Gaussian distribution
* :math:`x_{t+1}` is a linear transformation of :math:`x_t` and additive Gaussian noise
* :math:`z_{t+1}` is a linear transformation of :math:`x_{t+1}` and additive Gaussian noise

These assumptions imply that that :math:`x_t` is always a Gaussian distribution, even when :math:`z_t` is observed.  If this is the case, the distribution of `x_t|z_{1:T}` and `x_t | z_{1:T}` are completely specified by the covariance of the Gaussian distribution, namely its *mean* and *covariance*.  The Kalman Filter and Kalman Smoother calculate these values, respectively.

Formally, the Linear-Gaussian Model assumes that states and measurements are generated in the following way,

.. math:: x_0 \sim \text{Gaussian}(\mu_0, \Sigma_0)
.. math:: x_{t+1} = A x_t + b_t + \epsilon_{t}^{1}
.. math:: y_{t+1} = C x_{t+1} + d_t + \epsilon_{t}^2
.. math:: \epsilon_t^1 \sim \text{Gaussian}(0, Q)
.. math:: \epsilon_t^2 \sim \text{Gaussian}(0, R)

These assumptions mean that the Kalman Filter and Kalman Smoother work best if one is able to guess fairly well the vicinity of where the next state will be given the present, but cannot say *exactly* where it will be.  On the other hand, these methods will fail if there are multiple, disconnected areas where the next state could be, such as if a car turns at an intersection.

Filtering vs. Smoothing
=======================

Both Kalman Filtering and Kalman Smoothing aim to perform the same task: estimate the hidden state using the measurements.  Why then should one prefer one over the other?  The answer is that the Filter requires lower computational complexity while the Smoother gives better estimates.  Mathematically, the Filter only uses measurements :math:`z_1, \ldots, z_t` to estimates :math:`x_t`, but the Smoother uses :math:`z_1, \ldots, z_t, \ldots, z_T`.  In fact, the output of the Kalman Filter is necessary for implementing the Kalman Smoother.

In general, the computational complexity of the Kalman Filter and the Kalman Smoother is :math:`O(Td^3)`, and thus the Smoother should be preferred.  In practice, the Smoother takes roughly twice as the Filter as it must perform two passes over the measurements.  The only case where the Filter is better suited is when measurements :math:`z_t` come in a streaming fashion and estimates for `x_t` need to be updated online.

Usage
=====

The only class in this model is the :class:`KalmanFilter`.  It implements methods for Filtering, Smoothing, and EM all in one, as well as for sampling trajectories according to the generative model specified above.

The following snippet of code shows how one may sample from the model, learn model parameters, and apply the Kalman Filter and Kalman Smoother

.. literalinclude:: ../auto_examples/plot_kalman_simple.py
   :language: python

.. figure:: ../auto_examples/images/plot_kalman_simple_1.png
   :target: ../auto_examples/plot_kalman_simple.html
   :align: center
   :width: 400 px

.. topic:: Examples:

 * :ref:`example_plot_kalman.py`: Applying the Kalman Filter to the synthetic Kalman data set
