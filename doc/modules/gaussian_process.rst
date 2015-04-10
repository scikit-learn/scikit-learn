

.. _gaussian_process:

==================
Gaussian Processes
==================

.. currentmodule:: sklearn.gaussian_process

**Gaussian Processes for Machine Learning (GPML)** is a generic supervised
learning method primarily designed to solve *regression* problems. It has also
been extended to *probabilistic classification*, but in the present
implementation, this is only a post-processing of the *regression* exercise.

The advantages of Gaussian Processes for Machine Learning are:

    - The prediction interpolates the observations (at least for regular
      correlation models).

    - The prediction is probabilistic (Gaussian) so that one can compute
      empirical confidence intervals and exceedance probabilities that might be
      used to refit (online fitting, adaptive fitting) the prediction in some
      region of interest.

    - Versatile: different :ref:`linear regression models
      <linear_model>` and :ref:`correlation models
      <correlation_models>` can be specified. Common models are provided, but
      it is also possible to specify custom models provided they are
      stationary.

The disadvantages of Gaussian Processes for Machine Learning include:

    - It is not sparse. It uses the whole samples/features information to
      perform the prediction.

    - It loses efficiency in high dimensional spaces -- namely when the number
      of features exceeds a few dozens. It might indeed give poor performance
      and it loses computational efficiency.

    - Classification is only a post-processing, meaning that one first need
      to solve a regression problem by providing the complete scalar float
      precision output :math:`y` of the experiment one attempt to model.

Thanks to the Gaussian property of the prediction, it has been given varied
applications: e.g. for global optimization, probabilistic classification.


Gaussian Process Regression
===========================

Noise-free and noisy targets
----------------------------

.. figure:: ../auto_examples/gaussian_process/images/plot_gpr_noisy_targets_001.png
   :target: ../auto_examples/gaussian_process/plot_gpr_noisy_targets.html
   :align: center

.. figure:: ../auto_examples/gaussian_process/images/plot_gpr_noisy_targets_002.png
   :target: ../auto_examples/gaussian_process/plot_gpr_noisy_targets.html
   :align: center


Gaussian process regression (GPR) with noise-level estimation
-------------------------------------------------------------
This example illustrates that GPR with a sum-kernel including a WhiteKernel can
estimate the noise level of data. An illustration of the
log-marginal-likelihood (LML) landscape shows that there exist two local
maxima of LML.

.. figure:: ../auto_examples/gaussian_process/images/plot_gpr_noisy_000.png
   :target: ../auto_examples/gaussian_process/plot_gpr_noisy.html
   :align: center

The first corresponds to a model with a high noise level and a
large length scale, which explains all variations in the data by noise.

.. figure:: ../auto_examples/gaussian_process/images/plot_gpr_noisy_001.png
   :target: ../auto_examples/gaussian_process/plot_gpr_noisy.html
   :align: center

The second one has a smaller noise level and shorter length scale, which explains
most of the variation by the noise-free functional relationship. The second
model has a higher likelihood; however, depending on the initial value for the
hyperparameters, the gradient-based optimization might also converge to the
high-noise solution. It is thus important to repeat the optimization several
times for different initializations.

.. figure:: ../auto_examples/gaussian_process/images/plot_gpr_noisy_002.png
   :target: ../auto_examples/gaussian_process/plot_gpr_noisy.html
   :align: center


Comparison of GPR and Kernel Ridge Regression
---------------------------------------------

Both kernel ridge regression (KRR) and GPR learn
a target function by employing internally the "kernel trick". KRR learns a
linear function in the space induced by the respective kernel which corresponds
to a non-linear function in the original space. The linear function in the
kernel space is chosen based on the mean-squared error loss with
ridge regularization. GPR uses the kernel to define the covariance of
a prior distribution over the target functions and uses the observed training
data to define a likelihood function. Based on Bayes theorem, a (Gaussian)
posterior distribution over target functions is defined, whose mean is used
for prediction.

A major difference is that GPR can choose the kernel's hyperparameters based
on gradient-ascent on the marginal likelihood function while KRR needs to
perform a grid search on a cross-validated loss function (mean-squared error
loss). A further difference is that GPR learns a generative, probabilistic
model of the target function and can thus provide meaningful confidence
intervals and posterior samples along with the predictions while KRR only
provides predictions.

The following figure illustrates both methods on an artificial dataset, which
consists of a sinusoidal target function and strong noise. The figure compares
the learned model of KRR and GPR based on a ExpSineSquared kernel, which is
suited for learning periodic functions. The kernel's hyperparameters control
the smoothness (l) and periodicity of the kernel (p). Moreover, the noise level
of the data is learned explicitly by GPR by an additional WhiteKernel component
in the kernel and by the regularization parameter alpha of KRR.

.. figure:: ../auto_examples/gaussian_process/images/plot_compare_gpr_krr_001.png
   :target: ../auto_examples/gaussian_process/plot_compare_gpr_krr.html
   :align: center

The figure shows that both methods learn reasonable models of the target
function. GPR correctly identifies the periodicity of the function to be
roughly 2*pi (6.28), while KRR chooses the doubled periodicity 4*pi. Besides
that, GPR provides reasonable confidence bounds on the prediction which are not
available for KRR. A major difference between the two methods is the time
required for fitting and predicting: while fitting KRR is fast in principle,
the grid-search for hyperparameter optimization scales exponentially with the
number of hyperparameters ("curse of dimensionality"). The gradient-based
optimization of the parameters in GPR does not suffer from this exponential
scaling and is thus considerable faster on this example with 3-dimensional
hyperparameter space. The time for predicting is similar; however, generating
the variance of the predictive distribution of GPR takes considerable longer
than just predicting the mean.

Gaussian process regression (GPR) on Mauna Loa CO2 data
-------------------------------------------------------

This example is based on Section 5.4.3 of "Gaussian Processes for Machine
Learning" [RW2006]_. It illustrates an example of complex kernel engineering and
hyperparameter optimization using gradient ascent on the
log-marginal-likelihood. The data consists of the monthly average atmospheric
CO2 concentrations (in parts per million by volume (ppmv)) collected at the
Mauna Loa Observatory in Hawaii, between 1958 and 1997. The objective is to
model the CO2 concentration as a function of the time t.

The kernel is composed of several terms that are responsible for explaining
different properties of the signal:
 - a long term, smooth rising trend is to be explained by an RBF kernel. The
   RBF kernel with a large length-scale enforces this component to be smooth;
   it is not enforced that the trend is rising which leaves this choice to the
   GP. The specific length-scale and the amplitude are free hyperparameters.
 - a seasonal component, which is to be explained by the periodic
   ExpSineSquared kernel with a fixed periodicity of 1 year. The length-scale
   of this periodic component, controlling its smoothness, is a free parameter.
   In order to allow decaying away from exact periodicity, the product with an
   RBF kernel is taken. The length-scale of this RBF component controls the
   decay time and is a further free parameter.
 - smaller, medium term irregularities are to be explained by a
   RationalQuadratic kernel component, whose length-scale and alpha parameter,
   which determines the diffuseness of the length-scales, are to be determined.
   According to [RW2006], these irregularities can better be explained by
   a RationalQuadratic than an RBF kernel component, probably because it can
   accommodate several length-scales.
 - a "noise" term, consisting of an RBF kernel contribution, which shall
   explain the correlated noise components such as local weather phenomena,
   and a WhiteKernel contribution for the white noise. The relative amplitudes
   and the RBF's length scale are further free parameters.

Maximizing the log-marginal-likelihood after subtracting the target's mean
yields the following kernel with an LML of -84.483:

::

   2.5e+03 * RBF(l=49.8)
   + 6.68 * RBF(l=100) * ExpSineSquared(l=1.37, p=1)
   + 0.215 * RationalQuadratic(alpha=3.98, l=0.982)
   + 0.0381 * RBF(l=0.136) + WhiteKernel(c=0.0335)

Thus, most of the target signal (sqrt(2.5e+03)ppm = 50ppm) is explained by a
long-term rising trend (length-scale 49.8 years). The periodic component has
an amplitude of sqrt(6.68)ppm = 2.58ppm, a decay time of 100 years and a
length-scale of 1.37. The long decay time indicates that we have a locally very
close to periodic seasonal component. The correlated noise has an amplitude of
sqrt(0.0381)ppm = 0.195ppm with a length scale of 0.136 years and a white-noise
contribution of sqrt(0.0335)ppm = 0.183pm. Thus, the overall noise level is
very small, indicating that the data can be very well explained by the model.
The following figure shows also that the model makes very confident predictions
until around 2015.

.. figure:: ../auto_examples/gaussian_process/images/plot_gpr_co2_001.png
   :target: ../auto_examples/gaussian_process/plot_gpr_co2.html
   :align: center


Gaussian Process Classification
===============================

Probabilistic predictions with Gaussian process classification (GPC)
--------------------------------------------------------------------

This example illustrates the predicted probability of GPC for an RBF kernel
with different choices of the hyperparameters. The first figure shows the
predicted probability of GPC with arbitrarily chosen hyperparameters and with
the hyperparameters corresponding to the maximum log-marginal-likelihood (LML).

While the hyperparameters chosen by optimizing LML have a considerable larger
LML, they perform slightly worse according to the log-loss on test data. The
figure shows that this is because they exhibit a steep change of the class
probabilities at the class boundaries (which is good) but have predicted
probabilities close to 0.5 far away from the class boundaries (which is bad)
This undiesirable effect is caused by the Laplace approximation used
internally by GPC.

The second figure shows the log-marginal-likelihood for different choices of
the kernel's hyperparameters, highlighting the two choices of the
hyperparameters used in the first figure by black dots.

.. figure:: ../auto_examples/gaussian_process/images/plot_gpc_000.png
   :target: ../auto_examples/gaussian_process/plot_gpc.html
   :align: center

.. figure:: ../auto_examples/gaussian_process/images/plot_gpc_001.png
   :target: ../auto_examples/gaussian_process/plot_gpc.html
   :align: center


Iso-probability lines for Gaussian Processes classification (GPC)
-----------------------------------------------------------------

.. figure:: ../auto_examples/gaussian_process/images/plot_gpc_isoprobability_001.png
   :target: ../auto_examples/gaussian_process/plot_gpc_isoprobability.html
   :align: center


Illustration of Gaussian process classification (GPC) on the XOR dataset
------------------------------------------------------------------------

.. figure:: ../auto_examples/gaussian_process/images/plot_gpc_xor_001.png
   :target: ../auto_examples/gaussian_process/plot_gpc_xor.html
   :align: center


Kernels for Gaussian Processes
==============================
.. currentmodule:: sklearn.gaussian_process.kernels

Kernels (also called "covariance functions" in the context of GPs) are a crucial
ingredient of GPs which determine the shape of prior and posterior of the GP.
They encode the assumptions on the function being learned by defining the "similarity"
of two datapoints combined with the assumption that similar datapoints should
have similar target values. Two categories of kernels can be distinguished:
stationary kernels depend only on the distance of two datapoints and not on their
absolute values :math:`k(x_i, x_j)= k(d(x_i, x_j))` and are thus invariant to
translations in the input space, while non-stationary kernels
depend also on the specific values of the datapoints. Stationary kernels can further
be subdivided into isotropic and anisotropic kernels, where isotropic kernels are
also invariant to rotations in the input space. For more details, we refer to
Chapter 4 of [RW2006]_.

Basic kernels
-------------
The :class:`ConstantKernel` kernel can be used as part of a product-kernel
where it scales the magnitude of the other factor (kernel) or as part of a
sum-kernel, where it modifies the mean of the Gaussian process. It depends
on a parameter :math:`c`. It is defined as:

.. math::
   k(x_i, x_j) = c \;\forall\; x_1, x_2

The main use-case of the :class:`WhiteKernel` kernel is as part of a
sum-kernel where it explains the noise-component of the signal. Tuning its
parameter :math:`c` corresponds to estimating the noise-level.
It is defined as:

.. math::
    k(x_i, x_j) = c \text{ if } x_i == x_j \text{ else } 0


Kernel operators
----------------
Kernel operators take one or two base kernels and combine them into a new
kernel. The :class:`Sum` kernel takes two kernels :math:`k1` and :math:`k2`
and combines them via :math:`k_{sum}(X, Y) = k1(X, Y) + k2(X, Y)`.
The  :class:`Product` kernel takes two kernels :math:`k1` and :math:`k2`
and combines them via :math:`k_{product}(X, Y) = k1(X, Y) * k2(X, Y)`.
The :class:`Exponentiation` kernel takes one base kernel and a scalar parameter
:math:`exponent` and combines them via
:math:`k_{exp}(X, Y) = k(X, Y)^\text{exponent}`.

Radial-basis function (RBF) kernel
----------------------------------
The :class:`RBF` kernel is a stationary kernel. It is also known as the "squared
exponential" kernel. It is parameterized by a length-scale parameter :math:`l>0`, which
can either be a scalar (isotropic variant of the kernel) or a vector with the same
number of dimensions as the inputs :math:`x` (anisotropic variant of the kernel).
The kernel given by:

.. math::
   k(x_i, x_j) = \text{exp}\left(-\frac{1}{2} \vert d(x_i / l, x_j / l)^2\right)

This kernel is infinitely differentiable, which implies that GPs with this
kernel as covariance function have mean square derivatives of all orders, and are thus
very smooth. The prior and posterior of a GP resulting from an RBF kernel is shown in
the following figure:

.. figure:: ../auto_examples/gaussian_process/images/plot_gpr_prior_posterior_000.png
   :target: ../auto_examples/gaussian_process/plot_gpr_prior_posterior.html
   :align: center

Rational quadratic kernel
-------------------------

The :class:`RationalQuadratic` kernel can be seen as a scale mixture (an infinite sum)
of :class:`RBF` kernels with different characteristic length-scales. It is parameterized
by a length-scale parameter :math:`l>0` and a scale mixture parameter  :math:`\alpha>0`
Only the isotropic variant where :math:`l` is a scalar is supported at the moment.
The kernel given by:

.. math::
   k(x_i, x_j) = \left(1 + \frac{d(x_i, x_j)^2}{2\alpha l^2}\right)^\alpha

The prior and posterior of a GP resulting from an RBF kernel is shown in
the following figure:

.. figure:: ../auto_examples/gaussian_process/images/plot_gpr_prior_posterior_001.png
   :target: ../auto_examples/gaussian_process/plot_gpr_prior_posterior.html
   :align: center

Exp-Sine-Squared kernel
-----------------------

The :class:`ExpSineSquared` kernel allows modeling periodic functions.
It is parameterized by a length-scale parameter :math:`l>0` and a periodicity parameter
:math:`p>0`. Only the isotropic variant where :math:`l` is a scalar is supported at the moment.
The kernel given by:

.. math::
   k(x_i, x_j) = \text{exp}\left(-2 \text{sin}(\pi / p * d(x_i, x_j)) / l\right)^2

The prior and posterior of a GP resulting from an ExpSineSquared kernel is shown in
the following figure:

.. figure:: ../auto_examples/gaussian_process/images/plot_gpr_prior_posterior_002.png
   :target: ../auto_examples/gaussian_process/plot_gpr_prior_posterior.html
   :align: center

Dot-Product kernel
------------------

The :class:`DotProduct` kernel is non-stationary and can be obtained from linear regression
by putting :math:`N(0, 1)` priors on the coefficients of :math:`x_d (d = 1, . . . , D)` and
a prior of :math:`N(0, \sigma_0^2)` on the bias. The :class:`DotProduct` kernel is invariant to a rotation
of the coordinates about the origin, but not translations.
It is parameterized by a parameter :math:`\sigma_0^2`. For :math:`\sigma_0^2 = 0`, the kernel
is called the homogeneous linear kernel, otherwise it is inhomogeneous. The kernel is given by

.. math::
   k(x_i, x_j) = \sigma_0 ^ 2 + x_i \cdot x_j

The :class:`DotProduct` kernel is commonly combined with exponentiation. An example with exponent 2 is
shown in the following figure:

.. figure:: ../auto_examples/gaussian_process/images/plot_gpr_prior_posterior_003.png
   :target: ../auto_examples/gaussian_process/plot_gpr_prior_posterior.html
   :align: center


Legacy
======

An introductory regression example
----------------------------------

Say we want to surrogate the function :math:`g(x) = x \sin(x)`. To do so,
the function is evaluated onto a design of experiments. Then, we define a
GaussianProcess model whose regression and correlation models might be
specified using additional kwargs, and ask for the model to be fitted to the
data. Depending on the number of parameters provided at instantiation, the
fitting procedure may recourse to maximum likelihood estimation for the
parameters or alternatively it uses the given parameters.

.. figure:: ../auto_examples/gaussian_process/images/plot_gp_regression_001.png
   :target: ../auto_examples/gaussian_process/plot_gp_regression.html
   :align: center

::

    >>> import numpy as np
    >>> from sklearn import gaussian_process
    >>> def f(x):
    ...	    return x * np.sin(x)
    >>> X = np.atleast_2d([1., 3., 5., 6., 7., 8.]).T
    >>> y = f(X).ravel()
    >>> x = np.atleast_2d(np.linspace(0, 10, 1000)).T
    >>> gp = gaussian_process.GaussianProcess(theta0=1e-2, thetaL=1e-4, thetaU=1e-1)
    >>> gp.fit(X, y)  # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    GaussianProcess(beta0=None, corr=<function squared_exponential at 0x...>,
            normalize=True, nugget=array(2.22...-15),
            optimizer='fmin_cobyla', random_start=1, random_state=...
            regr=<function constant at 0x...>, storage_mode='full',
            theta0=array([[ 0.01]]), thetaL=array([[ 0.0001]]),
            thetaU=array([[ 0.1]]), verbose=False)
    >>> y_pred, sigma2_pred = gp.predict(x, eval_MSE=True)


Fitting Noisy Data
------------------

When the data to be fit includes noise, the Gaussian process model can be
used by specifying the variance of the noise for each point.
:class:`GaussianProcess` takes a parameter ``nugget`` which
is added to the diagonal of the correlation matrix between training points:
in general this is a type of Tikhonov regularization.  In the special case
of a squared-exponential correlation function, this normalization is
equivalent to specifying a fractional variance in the input.  That is

.. math::
   \mathrm{nugget}_i = \left[\frac{\sigma_i}{y_i}\right]^2

With ``nugget`` and ``corr`` properly set, Gaussian Processes can be
used to robustly recover an underlying function from noisy data:

.. figure:: ../auto_examples/gaussian_process/images/plot_gp_regression_002.png
   :target: ../auto_examples/gaussian_process/plot_gp_regression.html
   :align: center

.. topic:: Other examples

  * :ref:`example_gaussian_process_plot_gp_probabilistic_classification_after_regression.py`



Mathematical formulation
========================


The initial assumption
----------------------

Suppose one wants to model the output of a computer experiment, say a
mathematical function:

.. math::

        g: & \mathbb{R}^{n_{\rm features}} \rightarrow \mathbb{R} \\
           & X \mapsto y = g(X)

GPML starts with the assumption that this function is *a* conditional sample
path of *a* Gaussian process :math:`G` which is additionally assumed to read as
follows:

.. math::

        G(X) = f(X)^T \beta + Z(X)

where :math:`f(X)^T \beta` is a linear regression model and :math:`Z(X)` is a
zero-mean Gaussian process with a fully stationary covariance function:

.. math::

        C(X, X') = \sigma^2 R(|X - X'|)

:math:`\sigma^2` being its variance and :math:`R` being the correlation
function which solely depends on the absolute relative distance between each
sample, possibly featurewise (this is the stationarity assumption).

From this basic formulation, note that GPML is nothing but an extension of a
basic least squares linear regression problem:

.. math::

        g(X) \approx f(X)^T \beta

Except we additionally assume some spatial coherence (correlation) between the
samples dictated by the correlation function. Indeed, ordinary least squares
assumes the correlation model :math:`R(|X - X'|)` is one when :math:`X = X'`
and zero otherwise : a *dirac* correlation model -- sometimes referred to as a
*nugget* correlation model in the kriging literature.


The best linear unbiased prediction (BLUP)
------------------------------------------

We now derive the *best linear unbiased prediction* of the sample path
:math:`g` conditioned on the observations:

.. math::

    \hat{G}(X) = G(X | y_1 = g(X_1), ...,
                                y_{n_{\rm samples}} = g(X_{n_{\rm samples}}))

It is derived from its *given properties*:

- It is linear (a linear combination of the observations)

.. math::

    \hat{G}(X) \equiv a(X)^T y

- It is unbiased

.. math::

    \mathbb{E}[G(X) - \hat{G}(X)] = 0

- It is the best (in the Mean Squared Error sense)

.. math::

    \hat{G}(X)^* = \arg \min\limits_{\hat{G}(X)} \;
                                            \mathbb{E}[(G(X) - \hat{G}(X))^2]

So that the optimal weight vector :math:`a(X)` is solution of the following
equality constrained optimization problem:

.. math::

    a(X)^* = \arg \min\limits_{a(X)} & \; \mathbb{E}[(G(X) - a(X)^T y)^2] \\
                       {\rm s. t.} & \; \mathbb{E}[G(X) - a(X)^T y] = 0

Rewriting this constrained optimization problem in the form of a Lagrangian and
looking further for the first order optimality conditions to be satisfied, one
ends up with a closed form expression for the sought predictor -- see
references for the complete proof.

In the end, the BLUP is shown to be a Gaussian random variate with mean:

.. math::

    \mu_{\hat{Y}}(X) = f(X)^T\,\hat{\beta} + r(X)^T\,\gamma

and variance:

.. math::

    \sigma_{\hat{Y}}^2(X) = \sigma_{Y}^2\,
    ( 1
    - r(X)^T\,R^{-1}\,r(X)
    + u(X)^T\,(F^T\,R^{-1}\,F)^{-1}\,u(X)
    )

where we have introduced:

* the correlation matrix whose terms are defined wrt the autocorrelation
  function and its built-in parameters :math:`\theta`:

.. math::

    R_{i\,j} = R(|X_i - X_j|, \theta), \; i,\,j = 1, ..., m

* the vector of cross-correlations between the point where the prediction is
  made and the points in the DOE:

.. math::

    r_i = R(|X - X_i|, \theta), \; i = 1, ..., m

* the regression matrix (eg the Vandermonde matrix if :math:`f` is a polynomial
  basis):

.. math::

    F_{i\,j} = f_i(X_j), \; i = 1, ..., p, \, j = 1, ..., m

* the generalized least square regression weights:

.. math::

    \hat{\beta} =(F^T\,R^{-1}\,F)^{-1}\,F^T\,R^{-1}\,Y

* and the vectors:

.. math::

    \gamma & = R^{-1}(Y - F\,\hat{\beta}) \\
    u(X) & = F^T\,R^{-1}\,r(X) - f(X)

It is important to notice that the probabilistic response of a Gaussian Process
predictor is fully analytic and mostly relies on basic linear algebra
operations. More precisely the mean prediction is the sum of two simple linear
combinations (dot products), and the variance requires two matrix inversions,
but the correlation matrix can be decomposed only once using a Cholesky
decomposition algorithm.


The empirical best linear unbiased predictor (EBLUP)
----------------------------------------------------

Until now, both the autocorrelation and regression models were assumed given.
In practice however they are never known in advance so that one has to make
(motivated) empirical choices for these models :ref:`correlation_models`.

Provided these choices are made, one should estimate the remaining unknown
parameters involved in the BLUP. To do so, one uses the set of provided
observations in conjunction with some inference technique. The present
implementation, which is based on the DACE's Matlab toolbox uses the *maximum
likelihood estimation* technique -- see DACE manual in references for the
complete equations. This maximum likelihood estimation problem is turned into
a global optimization problem onto the autocorrelation parameters. In the
present implementation, this global optimization is solved by means of the
fmin_cobyla optimization function from scipy.optimize. In the case of
anisotropy however, we provide an implementation of Welch's componentwise
optimization algorithm -- see references.

For a more comprehensive description of the theoretical aspects of Gaussian
Processes for Machine Learning, please refer to the references below:

.. topic:: References:

    .. `DACE, A Matlab Kriging Toolbox
       <http://www2.imm.dtu.dk/~hbn/dace/>`_ S Lophaven, HB Nielsen, J
       Sondergaard 2002


    .. `Screening, predicting, and computer experiments
      <http://www.jstor.org/pss/1269548>`_ WJ Welch, RJ Buck, J Sacks,
      HP Wynn, TJ Mitchell, and MD Morris Technometrics 34(1) 15--25,
      1992


    .. [RW2006] `Gaussian Processes for Machine Learning
      <http://www.gaussianprocess.org/gpml/chapters/>`_ CE
      Rasmussen, CKI Williams, MIT Press, 2006 (Ed. T Diettrich)

    .. `The design and analysis of computer experiments
      <http://www.stat.osu.edu/~comp_exp/book.html>`_ TJ Santner, BJ
      Williams, W Notz Springer, 2003



.. _correlation_models:

Correlation Models
==================

Common correlation models matches some famous SVM's kernels because they are
mostly built on equivalent assumptions. They must fulfill Mercer's conditions
and should additionally remain stationary. Note however, that the choice of the
correlation model should be made in agreement with the known properties of the
original experiment from which the observations come. For instance:

* If the original experiment is known to be infinitely differentiable (smooth),
  then one should use the *squared-exponential correlation model*.
* If it's not, then one should rather use the *exponential correlation model*.
* Note also that there exists a correlation model that takes the degree of
  derivability as input: this is the Matern correlation model, but it's not
  implemented here (TODO).

For a more detailed discussion on the selection of appropriate correlation
models, see the book by Rasmussen & Williams in references.

.. _regression_models:


Regression Models
=================

Common linear regression models involve zero- (constant), first- and
second-order polynomials. But one may specify its own in the form of a Python
function that takes the features X as input and that returns a vector
containing the values of the functional set. The only constraint is that the
number of functions must not exceed the number of available observations so
that the underlying regression problem is not *underdetermined*.


Implementation details
======================

The present implementation is based on a translation of the DACE Matlab
toolbox.

.. topic:: References:

    * `DACE, A Matlab Kriging Toolbox
      <http://www2.imm.dtu.dk/~hbn/dace/>`_ S Lophaven, HB Nielsen, J
      Sondergaard 2002,

    * W.J. Welch, R.J. Buck, J. Sacks, H.P. Wynn, T.J. Mitchell, and M.D.
      Morris (1992). Screening, predicting, and computer experiments.
      Technometrics, 34(1) 15--25.
