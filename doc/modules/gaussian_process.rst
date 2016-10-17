

.. _gaussian_process:

==================
Gaussian Processes
==================

.. currentmodule:: sklearn.gaussian_process

**Gaussian Processes (GP)** are a generic supervised learning method designed
to solve *regression* and *probabilistic classification* problems.

The advantages of Gaussian processes are:

    - The prediction interpolates the observations (at least for regular
      kernels).

    - The prediction is probabilistic (Gaussian) so that one can compute
      empirical confidence intervals and decide based on those if one should
      refit (online fitting, adaptive fitting) the prediction in some
      region of interest.

    - Versatile: different :ref:`kernels
      <gp_kernels>` can be specified. Common kernels are provided, but
      it is also possible to specify custom kernels.

The disadvantages of Gaussian processes include:

    - They are not sparse, i.e., they use the whole samples/features information to
      perform the prediction.

    - They lose efficiency in high dimensional spaces -- namely when the number
      of features exceeds a few dozens.


.. _gpr:

Gaussian Process Regression (GPR)
=================================

.. currentmodule:: sklearn.gaussian_process

The :class:`GaussianProcessRegressor` implements Gaussian processes (GP) for
regression purposes. For this, the prior of the GP needs to be specified. The
prior mean is assumed to be constant and zero (for ``normalize_y=False``) or the
training data's mean (for ``normalize_y=True``). The prior's
covariance is specified by a passing a :ref:`kernel <gp_kernels>` object. The
hyperparameters of the kernel are optimized during fitting of
GaussianProcessRegressor by maximizing the log-marginal-likelihood (LML) based
on the passed ``optimizer``. As the LML may have multiple local optima, the
optimizer can be started repeatedly by specifying ``n_restarts_optimizer``. The
first run is always conducted starting from the initial hyperparameter values
of the kernel; subsequent runs are conducted from hyperparameter values
that have been chosen randomly from the range of allowed values.
If the initial hyperparameters should be kept fixed, `None` can be passed as
optimizer.

The noise level in the targets can be specified by passing it via the
parameter ``alpha``, either globally as a scalar or per datapoint.
Note that a moderate noise level can also be helpful for dealing with numeric
issues during fitting as it is effectively implemented as Tikhonov
regularization, i.e., by adding it to the diagonal of the kernel matrix. An
alternative to specifying the noise level explicitly is to include a
WhiteKernel component into the kernel, which can estimate the global noise
level from the data (see example below).

The implementation is based on Algorithm 2.1 of [RW2006]_. In addition to
the API of standard scikit-learn estimators, GaussianProcessRegressor:

* allows prediction without prior fitting (based on the GP prior)

* provides an additional method ``sample_y(X)``, which evaluates samples
  drawn from the GPR (prior or posterior) at given inputs

* exposes a method ``log_marginal_likelihood(theta)``, which can be used
  externally for other ways of selecting hyperparameters, e.g., via
  Markov chain Monte Carlo.


GPR examples
============

GPR with noise-level estimation
-------------------------------
This example illustrates that GPR with a sum-kernel including a WhiteKernel can
estimate the noise level of data. An illustration of the
log-marginal-likelihood (LML) landscape shows that there exist two local
maxima of LML.

.. figure:: ../auto_examples/gaussian_process/images/sphx_glr_plot_gpr_noisy_000.png
   :target: ../auto_examples/gaussian_process/plot_gpr_noisy.html
   :align: center

The first corresponds to a model with a high noise level and a
large length scale, which explains all variations in the data by noise.

.. figure:: ../auto_examples/gaussian_process/images/sphx_glr_plot_gpr_noisy_001.png
   :target: ../auto_examples/gaussian_process/plot_gpr_noisy.html
   :align: center

The second one has a smaller noise level and shorter length scale, which explains
most of the variation by the noise-free functional relationship. The second
model has a higher likelihood; however, depending on the initial value for the
hyperparameters, the gradient-based optimization might also converge to the
high-noise solution. It is thus important to repeat the optimization several
times for different initializations.

.. figure:: ../auto_examples/gaussian_process/images/sphx_glr_plot_gpr_noisy_002.png
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
the smoothness (length_scale) and periodicity of the kernel (periodicity).
Moreover, the noise level
of the data is learned explicitly by GPR by an additional WhiteKernel component
in the kernel and by the regularization parameter alpha of KRR.

.. figure:: ../auto_examples/gaussian_process/images/sphx_glr_plot_compare_gpr_krr_001.png
   :target: ../auto_examples/gaussian_process/plot_compare_gpr_krr.html
   :align: center

The figure shows that both methods learn reasonable models of the target
function. GPR correctly identifies the periodicity of the function to be
roughly :math:`2*\pi` (6.28), while KRR chooses the doubled periodicity
:math:`4*\pi` . Besides
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

GPR on Mauna Loa CO2 data
-------------------------

This example is based on Section 5.4.3 of [RW2006]_.
It illustrates an example of complex kernel engineering and
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
  According to [RW2006]_, these irregularities can better be explained by
  a RationalQuadratic than an RBF kernel component, probably because it can
  accommodate several length-scales.

- a "noise" term, consisting of an RBF kernel contribution, which shall
  explain the correlated noise components such as local weather phenomena,
  and a WhiteKernel contribution for the white noise. The relative amplitudes
  and the RBF's length scale are further free parameters.

Maximizing the log-marginal-likelihood after subtracting the target's mean
yields the following kernel with an LML of -83.214:

::

   34.4**2 * RBF(length_scale=41.8)
   + 3.27**2 * RBF(length_scale=180) * ExpSineSquared(length_scale=1.44,
                                                      periodicity=1)
   + 0.446**2 * RationalQuadratic(alpha=17.7, length_scale=0.957)
   + 0.197**2 * RBF(length_scale=0.138) + WhiteKernel(noise_level=0.0336)

Thus, most of the target signal (34.4ppm) is explained by a long-term rising
trend (length-scale 41.8 years). The periodic component has an amplitude of
3.27ppm, a decay time of 180 years and a length-scale of 1.44. The long decay
time indicates that we have a locally very close to periodic seasonal
component. The correlated noise has an amplitude of 0.197ppm with a length
scale of 0.138 years and a white-noise contribution of 0.197ppm. Thus, the
overall noise level is very small, indicating that the data can be very well
explained by the model. The figure shows also that the model makes very
confident predictions until around 2015

.. figure:: ../auto_examples/gaussian_process/images/sphx_glr_plot_gpr_co2_001.png
   :target: ../auto_examples/gaussian_process/plot_gpr_co2.html
   :align: center

.. _gpc:

Gaussian Process Classification (GPC)
=====================================

.. currentmodule:: sklearn.gaussian_process

The :class:`GaussianProcessClassifier` implements Gaussian processes (GP) for
classification purposes, more specifically for probabilistic classification,
where test predictions take the form of class probabilities.
GaussianProcessClassifier places a GP prior on a latent function :math:`f`,
which is then squashed through a link function to obtain the probabilistic
classification. The latent function :math:`f` is a so-called nuisance function,
whose values are not observed and are not relevant by themselves.
Its purpose is to allow a convenient formulation of the model, and :math:`f`
is removed (integrated out) during prediction. GaussianProcessClassifier
implements the logistic link function, for which the integral cannot be
computed analytically but is easily approximated in the binary case.

In contrast to the regression setting, the posterior of the latent function
:math:`f` is not Gaussian even for a GP prior since a Gaussian likelihood is
inappropriate for discrete class labels. Rather, a non-Gaussian likelihood
corresponding to the logistic link function (logit) is used.
GaussianProcessClassifier approximates the non-Gaussian posterior with a
Gaussian based on the Laplace approximation. More details can be found in
Chapter 3 of [RW2006]_.

The GP prior mean is assumed to be zero. The prior's
covariance is specified by a passing a :ref:`kernel <gp_kernels>` object. The
hyperparameters of the kernel are optimized during fitting of
GaussianProcessRegressor by maximizing the log-marginal-likelihood (LML) based
on the passed ``optimizer``. As the LML may have multiple local optima, the
optimizer can be started repeatedly by specifying ``n_restarts_optimizer``. The
first run is always conducted starting from the initial hyperparameter values
of the kernel; subsequent runs are conducted from hyperparameter values
that have been chosen randomly from the range of allowed values.
If the initial hyperparameters should be kept fixed, `None` can be passed as
optimizer.

:class:`GaussianProcessClassifier` supports multi-class classification
by performing either one-versus-rest or one-versus-one based training and
prediction.  In one-versus-rest, one binary Gaussian process classifier is
fitted for each class, which is trained to separate this class from the rest.
In "one_vs_one", one binary Gaussian process classifier is fitted for each pair
of classes, which is trained to separate these two classes. The predictions of
these binary predictors are combined into multi-class predictions. See the
section on :ref:`multi-class classification <multiclass>` for more details.

In the case of Gaussian process classification, "one_vs_one" might be
computationally  cheaper since it has to solve many problems involving only a
subset of the whole training set rather than fewer problems on the whole
dataset. Since Gaussian process classification scales cubically with the size
of the dataset, this might be considerably faster. However, note that
"one_vs_one" does not support predicting probability estimates but only plain
predictions. Moreover, note that :class:`GaussianProcessClassifier` does not
(yet) implement a true multi-class Laplace approximation internally, but
as discussed aboved is based on solving several binary classification tasks
internally, which are combined using one-versus-rest or one-versus-one.

GPC examples
============

Probabilistic predictions with GPC
----------------------------------

This example illustrates the predicted probability of GPC for an RBF kernel
with different choices of the hyperparameters. The first figure shows the
predicted probability of GPC with arbitrarily chosen hyperparameters and with
the hyperparameters corresponding to the maximum log-marginal-likelihood (LML).

While the hyperparameters chosen by optimizing LML have a considerable larger
LML, they perform slightly worse according to the log-loss on test data. The
figure shows that this is because they exhibit a steep change of the class
probabilities at the class boundaries (which is good) but have predicted
probabilities close to 0.5 far away from the class boundaries (which is bad)
This undesirable effect is caused by the Laplace approximation used
internally by GPC.

The second figure shows the log-marginal-likelihood for different choices of
the kernel's hyperparameters, highlighting the two choices of the
hyperparameters used in the first figure by black dots.

.. figure:: ../auto_examples/gaussian_process/images/sphx_glr_plot_gpc_000.png
   :target: ../auto_examples/gaussian_process/plot_gpc.html
   :align: center

.. figure:: ../auto_examples/gaussian_process/images/sphx_glr_plot_gpc_001.png
   :target: ../auto_examples/gaussian_process/plot_gpc.html
   :align: center


Illustration of GPC on the XOR dataset
--------------------------------------

.. currentmodule:: sklearn.gaussian_process.kernels

This example illustrates GPC on XOR data. Compared are a stationary, isotropic
kernel (:class:`RBF`) and a non-stationary kernel (:class:`DotProduct`). On this particular
dataset, the `DotProduct` kernel obtains considerably better results because the
class-boundaries are linear and coincide with the coordinate axes. In practice,
however, stationary kernels such as :class:`RBF` often obtain better results.

.. figure:: ../auto_examples/gaussian_process/images/sphx_glr_plot_gpc_xor_001.png
   :target: ../auto_examples/gaussian_process/plot_gpc_xor.html
   :align: center

.. currentmodule:: sklearn.gaussian_process


Gaussian process classification (GPC) on iris dataset
-----------------------------------------------------

This example illustrates the predicted probability of GPC for an isotropic
and anisotropic RBF kernel on a two-dimensional version for the iris-dataset.
This illustrates the applicability of GPC to non-binary classification.
The anisotropic RBF kernel obtains slightly higher log-marginal-likelihood by
assigning different length-scales to the two feature dimensions.

.. figure:: ../auto_examples/gaussian_process/images/sphx_glr_plot_gpc_iris_001.png
   :target: ../auto_examples/gaussian_process/plot_gpc_iris.html
   :align: center


.. _gp_kernels:

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

Gaussian Process Kernel API
---------------------------
The main usage of a :class:`Kernel` is to compute the GP's covariance between
datapoints. For this, the method ``__call__`` of the kernel can be called. This
method can either be used to compute the "auto-covariance" of all pairs of
datapoints in a 2d array X, or the "cross-covariance" of all combinations
of datapoints of a 2d array X with datapoints in a 2d array Y. The following
identity holds true for all kernels k (except for the :class:`WhiteKernel`):
``k(X) == K(X, Y=X)``

If only the diagonal of the auto-covariance is being used, the method ``diag()``
of a kernel can be called, which is more computationally efficient than the
equivalent call to ``__call__``: ``np.diag(k(X, X)) == k.diag(X)``

Kernels are parameterized by a vector :math:`\theta` of hyperparameters. These
hyperparameters can for instance control length-scales or periodicity of a
kernel (see below). All kernels support computing analytic gradients of
of the kernel's auto-covariance with respect to :math:`\theta` via setting
``eval_gradient=True`` in the ``__call__`` method. This gradient is used by the
Gaussian process (both regressor and classifier) in computing the gradient
of the log-marginal-likelihood, which in turn is used to determine the
value of :math:`\theta`, which maximizes the log-marginal-likelihood,  via
gradient ascent. For each hyperparameter, the initial value and the
bounds need to be specified when creating an instance of the kernel. The
current value of :math:`\theta` can be get and set via the property
``theta`` of the kernel object. Moreover, the bounds of the hyperparameters can be
accessed by the property ``bounds`` of the kernel. Note that both properties
(theta and bounds) return log-transformed values of the internally used values
since those are typically more amenable to gradient-based optimization.
The specification of each hyperparameter is stored in the form of an instance of
:class:`Hyperparameter` in the respective kernel. Note that a kernel using a
hyperparameter with name "x" must have the attributes self.x and self.x_bounds.

The abstract base class for all kernels is :class:`Kernel`. Kernel implements a
similar interface as :class:`Estimator`, providing the methods ``get_params()``,
``set_params()``, and ``clone()``. This allows setting kernel values also via
meta-estimators such as :class:`Pipeline` or :class:`GridSearch`. Note that due to the nested
structure of kernels (by applying kernel operators, see below), the names of
kernel parameters might become relatively complicated. In general, for a
binary kernel operator, parameters of the left operand are prefixed with ``k1__``
and parameters of the right operand with ``k2__``. An additional convenience
method is ``clone_with_theta(theta)``, which returns a cloned version of the
kernel but with the hyperparameters set to ``theta``. An illustrative example:

    >>> from sklearn.gaussian_process.kernels import ConstantKernel, RBF
    >>> kernel = ConstantKernel(constant_value=1.0, constant_value_bounds=(0.0, 10.0)) * RBF(length_scale=0.5, length_scale_bounds=(0.0, 10.0)) + RBF(length_scale=2.0, length_scale_bounds=(0.0, 10.0))
    >>> for hyperparameter in kernel.hyperparameters: print(hyperparameter)
    Hyperparameter(name='k1__k1__constant_value', value_type='numeric', bounds=array([[  0.,  10.]]), n_elements=1, fixed=False)
    Hyperparameter(name='k1__k2__length_scale', value_type='numeric', bounds=array([[  0.,  10.]]), n_elements=1, fixed=False)
    Hyperparameter(name='k2__length_scale', value_type='numeric', bounds=array([[  0.,  10.]]), n_elements=1, fixed=False)
    >>> params = kernel.get_params()
    >>> for key in sorted(params): print("%s : %s" % (key, params[key]))
    k1 : 1**2 * RBF(length_scale=0.5)
    k1__k1 : 1**2
    k1__k1__constant_value : 1.0
    k1__k1__constant_value_bounds : (0.0, 10.0)
    k1__k2 : RBF(length_scale=0.5)
    k1__k2__length_scale : 0.5
    k1__k2__length_scale_bounds : (0.0, 10.0)
    k2 : RBF(length_scale=2)
    k2__length_scale : 2.0
    k2__length_scale_bounds : (0.0, 10.0)
    >>> print(kernel.theta)  # Note: log-transformed
    [ 0.         -0.69314718  0.69314718]
    >>> print(kernel.bounds)  # Note: log-transformed
    [[       -inf  2.30258509]
     [       -inf  2.30258509]
     [       -inf  2.30258509]]


All Gaussian process kernels are interoperable with :mod:`sklearn.metrics.pairwise`
and vice versa: instances of subclasses of :class:`Kernel` can be passed as
``metric`` to pairwise_kernels`` from :mod:`sklearn.metrics.pairwise`. Moreover,
kernel functions from pairwise can be used as GP kernels by using the wrapper
class :class:`PairwiseKernel`. The only caveat is that the gradient of
the hyperparameters is not analytic but numeric and all those kernels support
only isotropic distances. The parameter ``gamma`` is considered to be a
hyperparameter and may be optimized. The other kernel parameters are set
directly at initialization and are kept fixed.


Basic kernels
-------------
The :class:`ConstantKernel` kernel can be used as part of a :class:`Product`
kernel where it scales the magnitude of the other factor (kernel) or as part
of a :class:`Sum` kernel, where it modifies the mean of the Gaussian process.
It depends on a parameter :math:`constant\_value`. It is defined as:

.. math::
   k(x_i, x_j) = constant\_value \;\forall\; x_1, x_2

The main use-case of the :class:`WhiteKernel` kernel is as part of a
sum-kernel where it explains the noise-component of the signal. Tuning its
parameter :math:`noise\_level` corresponds to estimating the noise-level.
It is defined as:e

.. math::
    k(x_i, x_j) = noise\_level \text{ if } x_i == x_j \text{ else } 0


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
The kernel is given by:

.. math::
   k(x_i, x_j) = \text{exp}\left(-\frac{1}{2} d(x_i / l, x_j / l)^2\right)

This kernel is infinitely differentiable, which implies that GPs with this
kernel as covariance function have mean square derivatives of all orders, and are thus
very smooth. The prior and posterior of a GP resulting from an RBF kernel are shown in
the following figure:

.. figure:: ../auto_examples/gaussian_process/images/sphx_glr_plot_gpr_prior_posterior_000.png
   :target: ../auto_examples/gaussian_process/plot_gpr_prior_posterior.html
   :align: center


Matérn kernel
-------------
The :class:`Matern` kernel is a stationary kernel and a generalization of the
:class:`RBF` kernel. It has an additional parameter :math:`\nu` which controls
the smoothness of the resulting function. It is parameterized by a length-scale parameter :math:`l>0`, which can either be a scalar (isotropic variant of the kernel) or a vector with the same number of dimensions as the inputs :math:`x` (anisotropic variant of the kernel). The kernel is given by:

.. math::

    k(x_i, x_j) = \sigma^2\frac{1}{\Gamma(\nu)2^{\nu-1}}\Bigg(\gamma\sqrt{2\nu} d(x_i / l, x_j / l)\Bigg)^\nu K_\nu\Bigg(\gamma\sqrt{2\nu} d(x_i / l, x_j / l)\Bigg),

As :math:`\nu\rightarrow\infty`, the Matérn kernel converges to the RBF kernel.
When :math:`\nu = 1/2`, the Matérn kernel becomes identical to the absolute
exponential kernel, i.e.,

.. math::
    k(x_i, x_j) = \sigma^2 \exp \Bigg(-\gamma d(x_i / l, x_j / l) \Bigg) \quad \quad \nu= \tfrac{1}{2}

In particular, :math:`\nu = 3/2`:

.. math::
    k(x_i, x_j) = \sigma^2 \Bigg(1 + \gamma \sqrt{3} d(x_i / l, x_j / l)\Bigg) \exp \Bigg(-\gamma \sqrt{3}d(x_i / l, x_j / l) \Bigg) \quad \quad \nu= \tfrac{3}{2}

and :math:`\nu = 5/2`:

.. math::
    k(x_i, x_j) = \sigma^2 \Bigg(1 + \gamma \sqrt{5}d(x_i / l, x_j / l) +\frac{5}{3} \gamma^2d(x_i / l, x_j / l)^2 \Bigg) \exp \Bigg(-\gamma \sqrt{5}d(x_i / l, x_j / l) \Bigg) \quad \quad \nu= \tfrac{5}{2}

are popular choices for learning functions that are not infinitely
differentiable (as assumed by the RBF kernel) but at least once (:math:`\nu =
3/2`) or twice differentiable (:math:`\nu = 5/2`).

The flexibility of controlling the smoothness of the learned function via :math:`\nu`
allows adapting to the properties of the true underlying functional relation.
The prior and posterior of a GP resulting from a Matérn kernel are shown in
the following figure:

.. figure:: ../auto_examples/gaussian_process/images/sphx_glr_plot_gpr_prior_posterior_004.png
   :target: ../auto_examples/gaussian_process/plot_gpr_prior_posterior.html
   :align: center

See [RW2006]_, pp84 for further details regarding the
different variants of the Matérn kernel.

Rational quadratic kernel
-------------------------

The :class:`RationalQuadratic` kernel can be seen as a scale mixture (an infinite sum)
of :class:`RBF` kernels with different characteristic length-scales. It is parameterized
by a length-scale parameter :math:`l>0` and a scale mixture parameter  :math:`\alpha>0`
Only the isotropic variant where :math:`l` is a scalar is supported at the moment.
The kernel is given by:

.. math::
   k(x_i, x_j) = \left(1 + \frac{d(x_i, x_j)^2}{2\alpha l^2}\right)^\alpha

The prior and posterior of a GP resulting from an RBF kernel are shown in
the following figure:

.. figure:: ../auto_examples/gaussian_process/images/sphx_glr_plot_gpr_prior_posterior_001.png
   :target: ../auto_examples/gaussian_process/plot_gpr_prior_posterior.html
   :align: center

Exp-Sine-Squared kernel
-----------------------

The :class:`ExpSineSquared` kernel allows modeling periodic functions.
It is parameterized by a length-scale parameter :math:`l>0` and a periodicity parameter
:math:`p>0`. Only the isotropic variant where :math:`l` is a scalar is supported at the moment.
The kernel is given by:

.. math::
   k(x_i, x_j) = \text{exp}\left(-2 \text{sin}(\pi / p * d(x_i, x_j)) / l\right)^2

The prior and posterior of a GP resulting from an ExpSineSquared kernel are shown in
the following figure:

.. figure:: ../auto_examples/gaussian_process/images/sphx_glr_plot_gpr_prior_posterior_002.png
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

.. figure:: ../auto_examples/gaussian_process/images/sphx_glr_plot_gpr_prior_posterior_003.png
   :target: ../auto_examples/gaussian_process/plot_gpr_prior_posterior.html
   :align: center

References
----------

    * `[RW2006]
      <http://www.gaussianprocess.org/gpml/chapters/>`_
      **Gaussian Processes for Machine Learning**,
      Carl Eduard Rasmussen and Christopher K.I. Williams, MIT Press 2006.
      Link to an official complete PDF version of the book
      `here <http://www.gaussianprocess.org/gpml/chapters/RW.pdf>`_ .

.. currentmodule:: sklearn.gaussian_process




Legacy Gaussian Processes
=========================

In this section, the implementation of Gaussian processes used in scikit-learn
until release 0.16.1 is described. Note that this implementation is deprecated
and will be removed in version 0.18.

An introductory regression example
----------------------------------

Say we want to surrogate the function :math:`g(x) = x \sin(x)`. To do so,
the function is evaluated onto a design of experiments. Then, we define a
GaussianProcess model whose regression and correlation models might be
specified using additional kwargs, and ask for the model to be fitted to the
data. Depending on the number of parameters provided at instantiation, the
fitting procedure may recourse to maximum likelihood estimation for the
parameters or alternatively it uses the given parameters.


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
used to robustly recover an underlying function from noisy data.



Mathematical formulation
------------------------


The initial assumption
^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

.. _correlation_models:

Correlation Models
------------------

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
-----------------

Common linear regression models involve zero- (constant), first- and
second-order polynomials. But one may specify its own in the form of a Python
function that takes the features X as input and that returns a vector
containing the values of the functional set. The only constraint is that the
number of functions must not exceed the number of available observations so
that the underlying regression problem is not *underdetermined*.


Implementation details
----------------------

The implementation is based on a translation of the DACE Matlab
toolbox.

.. topic:: References:

    * `DACE, A Matlab Kriging Toolbox
      <http://imedea.uib-csic.es/master/cambioglobal/Modulo_V_cod101615/Lab/lab_maps/krigging/DACE-krigingsoft/dace/dace.pdf>`_ S Lophaven, HB Nielsen, J
      Sondergaard 2002,

    * W.J. Welch, R.J. Buck, J. Sacks, H.P. Wynn, T.J. Mitchell, and M.D.
      Morris (1992). Screening, predicting, and computer experiments.
      Technometrics, 34(1) 15--25.
