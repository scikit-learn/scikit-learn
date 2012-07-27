.. _kalman:

=============
Kalman Filter
=============

.. currentmodule:: sklearn.kalman

The **Kalman Filter** is the de facto standard algorithm for tracking a single
moving object in discrete time.  Given a sequence of position observations, the
Kalman Filter is automatically able to infer model parameters as well as other
linearly-related variables such as velocity and acceleration.

.. figure:: ../auto_examples/kalman/images/plot_sin_1.png
   :target: ../auto_examples/kalman/plot_sin.html
   :width: 600 px
   :align: center

Like the Hidden Markov Model, the Kalman Filter estimates a sequence of hidden
states by performing Bayesian Inference on a *generative probabilistic model*
given a sequence of measurements.  Unlike the Hidden Markov Model, however, the
Kalman Filter is designed to work with continuous state and observation spaces.

The advantages of Kalman Filter are:

    + No need to provide labeled training data

    + Ability to handle noisy observations

The disadvantages are:

    + Computational complexity is cubic in the size of the state space

    + Parameter optimization is non-convex and can thus only find local optima.

    + Inability to cope with non-Gaussian noise

Usage
=====

While :class:`KalmanFilter` is the only class in this module, there are
actually three different algorithms implemented.  The first two are the Kalman
Filter and Kalman Smoother, two algorithms for estimating the hidden target
states.  The third is the EM algorithm, which makes use the both of the above
to estimate the model's parameters.

In order to make use of these algorithms, one need only supply an initial guess
for a subset of the model's parameters, one defining the size of the state
space and another the size of the measurement space::

    >>> from sklearn.kalman import KalmanFilter
    >>> import numpy as np
    >>> kf = KalmanFilter(transition_covariance=0.1, observation_covariance=np.eye(2))


We may then fit the Kalman Filter to a sequence of observations, and use it to
predict the underlying hidden states of any observation sequence::

    >>> measurements = [[1,0], [0,0], [0,1]]
    >>> kf.fit(measurements).predict([[2,0], [2,1], [2,2]])
    array([[ 0.75617284],
           [ 1.26661721],
           [ 1.53838218]])

The Kalman Filter is parameterized by 3 arrays for state transitions, 3 for
measurements, and 2 more for initial conditions.  Their names and function are
described in the next section.

.. topic:: Examples:

 * :ref:`example_kalman_plot_sin.py`


Mathematical Formulation
========================

In order to understand when the algorithms in this module will be effective, it
is important to understand what assumptions are being made.  To make notation
concise,  we refer to the hidden states as :math:`x_t`, the measurements as
:math:`z_t`, and the parameters of the :class:`KalmanFilter` class as follows,

    +----------------------------+------------------+
    |    Parameter Name          |      Notation    |
    +----------------------------+------------------+
    | `initial_state_mean`       | :math:`\mu_0`    |
    +----------------------------+------------------+
    | `initial_state_covariance` | :math:`\Sigma_0` |
    +----------------------------+------------------+
    | `transition_matrices`      | :math:`A`        |
    +----------------------------+------------------+
    | `transition_offsets`       | :math:`b`        |
    +----------------------------+------------------+
    | `transition_covariance`    | :math:`Q`        |
    +----------------------------+------------------+
    | `observation_matrices`     | :math:`C`        |
    +----------------------------+------------------+
    | `observation_offsets`      | :math:`b`        |
    +----------------------------+------------------+
    | `observation_covariance`   | :math:`R`        |
    +----------------------------+------------------+

In words, the Linear-Gaussian model assumes that for all time steps :math:`t =
0, \ldots, T-1` (here, :math:`T` is the number of time steps),

* :math:`x_0` is distributed according to a Gaussian distribution
* :math:`x_{t+1}` is an affine transformation of :math:`x_t` and additive
  Gaussian noise
* :math:`z_{t}` is an affine transformation of :math:`x_{t}` and additive
  Gaussian noise

.. figure:: ../auto_examples/kalman/images/plot_pomp_diagram_1.png
   :target: ../auto_examples/kalman/plot_pomp_diagram.html
   :width: 350 px
   :align: center

These assumptions imply that that :math:`x_t` is always a Gaussian
distribution, even when :math:`z_t` is observed.  If this is the case, the
distribution of :math:`x_t|z_{1:t}` and :math:`x_t | z_{1:T-1}` are completely
specified by the parameters of the Gaussian distribution, namely its *mean* and
*covariance*.  The Kalman Filter and Kalman Smoother calculate these values,
respectively.

Formally, the Linear-Gaussian Model assumes that states and measurements are
generated in the following way,

.. math::

    x_0               & \sim \text{Gaussian}(\mu_0, \Sigma_0)   \\
    x_{t+1}           & = A_t x_t + b_t + \epsilon_{t+1}^{1}    \\
    y_{t}             & = C_t x_t + d_t + \epsilon_{t}^2        \\
    \epsilon_t^1      & \sim \text{Gaussian}(0, Q)              \\
    \epsilon_{t}^2    & \sim \text{Gaussian}(0, R)

The Gaussian distribution is characterized by its single mode and exponentially
decreasing tails, meaning that the Kalman Filter and Kalman Smoother work best
if one is able to guess fairly well the vicinity of the next state given the
present, but cannot say *exactly* where it will be.  On the other hand, these
methods will fail if there are multiple, disconnected areas where the next
state could be, such as if a car turns one of three ways at an intersection.

.. topic:: References:

 * Abbeel, Pieter. "Maximum Likelihood, EM".
   http://www.cs.berkeley.edu/~pabbeel/cs287-fa11/
 * Yu, Byron M. and Shenoy, Krishna V. and Sahani, Maneesh. "Derivation of
   Kalman Filtering and Smoothing Equations".
   http://www.ece.cmu.edu/~byronyu/papers/derive_ks.pdf
 * Ghahramani, Zoubin and Hinton, Geoffrey E. "Parameter Estimation for
   Linear Dynamical Systems."
   http://mlg.eng.cam.ac.uk/zoubin/course04/tr-96-2.pdf
 * Welling, Max. "The Kalman Filter".
   http://www.cs.toronto.edu/~welling/classnotes/papers_class/KF.ps.gz

Filtering vs. Smoothing
=======================

Both Kalman Filtering and Kalman Smoothing aim to perform the same task:
estimate the hidden states using the measurements.  Why then should one prefer
one over the other?  The answer is that the Filter requires lower computational
complexity while the Smoother gives better estimates.  Mathematically, the
Filter calculates the mean and covariance of the Normal distribution
representing :math:`P(x_t | z_0, \ldots, z_t)` while the Smoother calculates
the same for :math:`P(x_t | z_0, \ldots, z_T)`; that is, all measurements.  In
fact, the output of the Kalman Filter is necessary to implement the Kalman
Smoother.

In general, the computational complexity of the Kalman Filter and the Kalman
Smoother are both :math:`O(Td^3)` where :math:`T` is the total number of time
steps and :math:`d` is the dimensionality of the state space, and thus the
Smoother should be preferred.  In practice, the Smoother takes roughly twice as
long as the Filter as it must perform two passes over the measurements.  The
only case where the Filter is better suited is when measurements :math:`z_t`
come in a streaming fashion and estimates for :math:`x_t` need to be updated
online.  If that is the case, one should look at
:func:`KalmanFilter.filter_update`.

Finally, textbook examples of the Kalman Filter and Kalman Smoother often
assume :math:`x_t` ranges from :math:`t = 0 \ldots T` while :math:`z_t` ranges
from :math:`t = 1 \ldots T`.  This module assumes both :math:`x_t` and
:math:`z_t` range from :math:`t = 0 \ldots T-1`.

.. topic:: Examples:

 * :ref:`example_kalman_plot_online.py`
 * :ref:`example_kalman_plot_filter.py`


EM Algorithm
============

The Expectation-Maximization Algorithm, better known as the EM algorithm, is
actually a bit of a misnomer; it is more like an algorithm *template* than an
algorithm in and of itself.  The EM algorithm seeks to maximize the likelihood
of all measurements by estimating the unknown variables (in this case, the
hidden states) with a fixed set of parameters, then maximizing the expected
value of the log likelihood of all measurements with respect to the parameters,
and repeating.  In mathematical notation, if we define :math:`\theta = (A, b,
C, d, Q, R, \mu_0, \Sigma_0)`, then the EM Algorithm works by finding a closed
form expression for,

.. math::

    P(x_{0:T-1} | \theta_{\text{old}}, z_{0:T-1})

then uses that expression to find,

.. math::

    \theta_{\text{new}} = \arg\max_{\theta'} \mathbb{E}_{x_{0:T-1}} [
      \log P(x_{0:T-1}, z_{0:T-1} | \theta')  | z_{0:T-1}, \theta_{\text{old}}
    ]

:math:`\theta_{\text{new}}` then takes the place of :math:`\theta_{\text{old}}`
and the process is repeated.  On a practical note, each iteration of the EM
algorithm requires running the Kalman Smoother, and thus the running time is
:math:`O(n_{\text{iter}} T d^3)` when the EM algorithm is run for
:math:`n_{\text{iter}}` iterations.

The EM algorithm is always guaranteed to converge, but not necessarily to a
global optimum.  Thus, it is important to start with a good guess for the
original parameter values (typically an order of magnitude is sufficient).
Secondly, the EM algorithm implemented here does not support regularization, so
values parameters can grow extremely out of hand with insufficient data.

.. topic:: Examples:

 * :ref:`example_kalman_plot_em.py`


Missing Observations
====================

A real system will often get measurements at regular points in time, but there
will also be times when the sensor fails.  :mod:`sklearn.kalman` offers you the
ability to continue applying all of its implemented algorithms even if this is
the case.  In order to use this feature, one simply needs to wrap the
measurements in :mod:`numpy.ma` and mark a timestep as masked::

  >>> from sklearn.datasets import load_kalman_data
  >>> import numpy as np
  >>> import numpy.ma as ma
  >>> Z = load_kalman_data().data
  >>> Z = ma.array(Z, mask=np.zeros(Z.shape))
  >>> Z[5] = ma.masked  # observation at time step 5 will now be ignored
  >>> Z[5] = ma.nomask  # observation will be recognized again

.. topic:: Examples:

 * :ref:`example_kalman_plot_missing.py`
