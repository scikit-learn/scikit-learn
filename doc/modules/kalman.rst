.. _kalman:

=============
Kalman Filter
=============

.. currentmodule:: sklearn.kalman

The Kalman Filter is a unsupervised algorithm for tracking a single object in a
continuous state space.  Given a sequence of noisy measurements, the Kalman
Filter is able to recover the "true state" of the underling object being
tracked. Common uses for the Kalman Filter include radar and sonar tracking and
state estimation in robotics.

The advantages of Kalman Filter are:

    + No need to provide labeled training data

    + Ability to handle noisy observations

The disadvantages are:

    + Computational complexity is cubic in the size of the state space

    + Parameter optimization is non-convex and can thus only find local optima

    + Inability to cope with non-Gaussian noise

.. figure:: ../auto_examples/kalman/images/plot_sin_1.png
   :target: ../auto_examples/kalman/plot_sin.html
   :width: 600 px
   :align: center


Usage
=====

This module implements two algorithms for tracking: the Kalman Filter and
Kalman Smoother.  In addition, model parameters which are traditionally
specified by hand can also be learned by the implemented EM algorithm without
any labeled training data.  All three algorithms are contained in the
:class:`KalmanFilter` class in this module.

In order to apply the Kalman Smoother, one need only specify the size of the
state and observation space.  This can be done directly by setting
:attr:`n_dim_state` or :attr:`n_dim_obs` or indirectly by specifying an initial
value for any of the model parameters from which the former can be derived::

    >>> from sklearn.kalman import KalmanFilter
    >>> import numpy as np
    >>> kf = KalmanFilter(initial_state_mean=0, n_dim_obs=2)

The traditional Kalman Filter assumes that model parameters are known
beforehand.  The :class:`KalmanFilter` class however can learn parameters using
:func:`KalmanFilter.fit` (fitting is optional).  Then the hidden sequence of
states can be predicted using :func:`KalmanFilter.predict`::

    >>> measurements = [[1,0], [0,0], [0,1]]
    >>> kf.fit(measurements).predict([[2,0], [2,1], [2,2]])
    array([[ 0.85819709],
           [ 1.77811829],
           [ 2.19537816]])

The Kalman Filter is parameterized by 3 arrays for state transitions, 3 for
measurements, and 2 more for initial conditions.  Their names and function are
described in the next section.

.. topic:: Examples:

 * :ref:`example_kalman_plot_sin.py`


Parameter Selection
-------------------

Unlike most other algorithms, the Kalman Filter and Kalman Smoother are
traditionally used with parameters already given. The :class:`KalmanFilter`
class can thus be initialized with any subset of the usual model parameters and
used without fitting. Sensible defaults values are given for all unspecified
parameters (zeros for all 1-dimensional arrays and identity matrices for all
2-dimensional arrays).

A Kalman Filter/Smoother is fully specified by its initial conditions
(:attr:`initial_state_mean` and :attr:`initial_state_covariance`), its
transition parameters (:attr:`transition_matrices`, :attr:`transition_offsets`,
:attr:`transition_covariance`), and its observation parameters
(:attr:`observation_matrices`, :attr:`observation_offsets`,
:attr:`observation_covariance`). These parameters define a probabilistic model
from which the unobserved states and observed measurements are assumed to be
sampled from. The following code illustrates in one dimension what this process
is.

.. code-block:: python

    from scipy.states import norm
    import numpy as np
    states = np.zeros((n_timesteps, n_dim_state))
    measurements = np.zeros((n_timesteps, n_dim_obs))
    for t in range(n_timesteps-1):
       if t == 0:
          states[t] = norm.rvs(initial_state_mean, np.sqrt(initial_state_covariance))
          measurements[t] = (
              np.dot(observation_matrices, states[t])
              + observation_offsets
              + norm.rvs(0, np.sqrt(observation_covariance))
          )
      states[t+1] = (
          np.dot(transition_matrices, states[t])
          + transition_offsets
          + norm.rvs(0, np.sqrt(transition_covariance))
      )
      measurements[t+1] = (
          np.dot(observation_matrices, states[t+1])
          + observation_offsets
          + norm.rvs(np.sqrt(observation_covariance))
      )

The selection of these variables is not an easy one, and, as shall be explained
in the section on fitting, should not be left to :func:`KalmanFilter.fit`
alone. If one ignores the random noise, the parameters dictate that *the next
state and the current measurement should be an affine function of the current
state*. The additive noise term is then simply a way to deal with unaccounted
error.

A simple example to illustrate the model parameters is a free falling ball in
one dimension. The state vector can be represented by the position, velocity,
and acceleration of the ball, and the transition matrix is defined by the
equation::

    position[t+dt] = position[t] + velocity[t] dt + 0.5 acceleration[t] dt^2

Taking the zeroth, first, and second derivative of the above equation with
respect to `dt` gives the rows of transition matrix. We may also set the
transition offset to zero for the position and velocity components and -9.8
for the acceleration component in order to account for gravity's pull.

It is often very difficult to guess what appropriate values are for for the
transition and observation covariance, so it is common to use some constant
multiplied by the identity matrix. Increasing this constant is equivalent to
saying you believe there is more noise in the system. This constant is the
amount of variance you expect to see along each dimensiona during state
transitions and measurements, respectively.


Prediction
----------

The :class:`KalmanFilter` class comes equipped with two algorithms for
prediction: the Kalman Filter and the Kalman Smoother. While the former can be
updated recursively (making it ideal for online state estimation), the latter
can only be done in batch. These two algorithms are accessible via
:func:`KalmanFilter.filter`, :func:`KalmanFilter.filter_update`, and
:func:`KalmanFilter.predict`.

Functionally, Kalman Smoother should always be preferred. Unlike the Kalman
Filter, the Smoother is able to incorporate "future" measurements as well as
past ones at the same computational cost of :math:`O(Td^3)` where :math:`T` is
the number of time steps and `d` is the dimensionality of the state space. The
only reason to prefer the Kalman Filter over the Smoother is in its ability to
incorporate new measurements in an online manner::

    means, covariances = kf.filter(measurements)
    next_mean, next_covariance = kf.filter_update(
        means[-1], covariances[-1], new_measurement
    )

Both the Kalman Filter and Kalman Smoother are able to use parameters which
vary with time.  In order to use this, one need only pass in an array
:attr:`n_timesteps` in length along its first axis::

    >>> transition_offsets = [[-1], [0], [1], [2]]
    >>> kf = KalmanFilter(transition_offsets=transition_offsets, n_dim_obs=1)

.. topic:: Examples:

 * :ref:`example_kalman_plot_online.py`
 * :ref:`example_kalman_plot_filter.py`


Fitting
-------

In addition to the Kalman Filter and Kalman Smoother, the :class:`KalmanFilter`
class implements the Expectation-Maximization algorithm. This iterative
algorithm is a way to maximize the likelihood of the observed measurements
(recall the probabilistic model induced by the model parameters), which is
unfortunately a non-convex optimization problem. This means that even when the
EM algorithm converges, there is no guarantee that it has converged to an
optimal value. Thus it is important to select good initial parameter values.

A second consideration when using the EM algorithm is that the algorithm lacks
regularization, meaning that parameter values may diverge to infinity in order
to make the measurements more likely. Thus it is important to choose *which*
parameters to optimize via the :attr:`em_vars` parameter of
:class:`KalmanFilter`.  For example, in order to only optimize the transition
and observation covariance matrices, one may instantiate :class:`KalmanFilter`
like so::

    >>> kf = KalmanFilter(em_vars=['transition_covariance', 'observation_covariance'])

It is customary optimize only the :attr:`transition_covariance`,
:attr:`observation_covariance`, :attr:`initial_state_mean`, and
:attr:`initial_state_covariance`, which is the default when :attr:`em_vars` is
unspecified. In order to avoid overfitting, it is also possible to specify the
number of iterations of the EM algorithm to run during fitting::

    >>> kf.fit(X, n_iter=5)

Each iteration of the EM algorithm requires running the Kalman Smoother anew,
so its computational complexity is :math:`O(Tnd^3)` where :math:`T` is the
number of time steps, `n` is the number of iterations, and `d` is the size of
the state space.

.. topic:: Examples:

 * :ref:`example_kalman_plot_em.py`


Missing Measurements
--------------------

In real world systems, it is common to have sensors occasionally fail.  The
Kalman Filter, Kalman Smoother, and EM algorithm are all equipped to handle
this scenario. To make use of it, one only need apply a NumPy mask to the
measurement at the missing time step::

    >>> from numpy import ma
    >>> X = ma.array([1,2,3])
    >>> X[1] = ma.masked  # hide measurement at time step 1
    >>> kf.fit(X).predict(X)

.. topic:: Examples:

 * :ref:`example_kalman_plot_missing.py`


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
