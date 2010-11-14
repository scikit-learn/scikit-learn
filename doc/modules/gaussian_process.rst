=======================================
Gaussian Processes for Machine Learning
=======================================

.. currentmodule:: scikits.learn.gpml

**Gaussian Processes for Machine Learning (GPML)** is a supervised learning
method used for *regression*. It can also be used for *probabilistic classification*, 
but it is only a post-processing of the *regression* exercise.

The advantages of Gaussian Processes for Machine Learning are:
      
    - The prediction interpolates the observations (at least for regular
      correlation models).

    - The prediction is probabilistic (Gaussian) so that one can compute 
      empirical confidence intervals and exceedence probabilities that might be 
      used to refit (online fitting, adaptive fitting) the prediction in some 
      region of interest.

    - Versatile: different :ref:`linear regression models <linear_regression_models>` and
      :ref:`correlation models <correlation_models>` can be specified. Common models are
      provided, but it is also possible to specify custom models.

The disadvantages of Gaussian Processes for Machine Learning include:

    - It is not sparse. It uses the whole samples/features information to 
      perform the prediction.

    - It loses efficiency in high dimensional spaces -- namely when the number
      of features exceeds a few dozens. It might indeed give poor performance
      and it becomes computationally inefficient.

    - Classification is only a post-processing, meaning that one first need
      to solve a regression problem by providing the complete scalar float
      precision output :math:`y` of the computer experiment one attempt to model.

Thanks to the Gaussian property of the prediction, it has been given varied 
applications: e.g. for global optimization, probabilistic classification.


Mathematical formulation
========================

The initial assumption
~~~~~~~~~~~~~~~~~~~~~~

Suppose one wants to model the output of a computer experiment, say a
mathematical function:

.. math::

        g: & \mathbb{R}^{n_{\rm features}} \rightarrow \mathbb{R} \\
           & X \mapsto y = g(X)

GPML starts with the assumption that this function is a conditionnal sample path
of a Gaussian process :math:`G` which is additionally assumed to read as follows:

.. math::

        G(X) = f(X)^T \beta + Z(X)

where :math:`f(X)^T \beta` is a linear regression model and :math:`Z(X)` is a zero-mean Gaussian
process with a fully stationary covariance function:

.. math::

        C(X, X') = \sigma^2 R(|X - X'|)

:math:`\sigma^2` being its variance and :math:`R` being the correlation function which solely
depends on the absolute relative distance between each sample -- possibly featurewise.

From this basic formulation, note that GPML is nothing but an extension of a
basic least squares linear regression problem:

.. math::

        g(X) \approx f(X)^T \beta

Except we additionaly assume some spatial coherence (correlation) between the 
samples dictated by the correlation function. Indeed, ordinary least squares assumes the
correlation model :math:`R(|X - X'|)` is one when :math:`X = X'` and zero otherwise : a *dirac* correlation
model -- sometimes referred to as a *nugget* correlation model in the kriging literature.


The best linear unbiased prediction (BLUP)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We now derive the *best linear unbiased prediction* of the sample path :math:`g`
conditioned by the observations:

.. math::

    \hat{G}(X) \sim G(X | y_1 = g(X_1), ..., y_{n_{\rm samples}} = g(X_{n_{\rm samples}}))

It is derived from its *given properties*:

- It is linear (a linear combination of the observations)

.. math::

    \hat{G}(X) \equiv a(X)^T y

- It is unbiased

.. math::

    \mathbb{E}[G(X) - \hat{G}(X)] = 0

- It is the best (in the Mean Squared Error sense)

.. math::

    \hat{G}(X)^* = \arg \min\limits_{\hat{G}(X)} \; \mathbb{E}[(G(X) - \hat{G}(X))^2]

So that the optimal weight vector :math:`a(X)` is solution of the following equality constrained optimization
problem:

.. math::

    a(X)^* = \arg \min\limits_{a(X)} & \; \mathbb{E}[(G(X) - a(X)^T y)^2] \\
                       {\rm s. t.} & \; \mathbb{E}[G(X) - a(X)^T y] = 0

Rewriting this constrained optimization problem in the form of a Lagrangian and looking further for
the first order optimality conditions to be satisfied, one ends up with a closed form expression for the 
sought predictor -- see references for the complete proof.

In the end, the BLUP is shown to be a Gaussian random variate whose moments expressions are given in reference.
   


The empirical best linear unbiased predictor (EBLUP)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Until now, both the autocorrelation and regression models were assumed given. In 
practice however they are never known in advance so that one has to make (motivated)
empirical choices for these models :ref:`correlation_models`.

Provided these choices are made, one should estimate the remaining unknown 
parameters involved in the BLUP. To do so, one uses the set of provided observations 
in conjunction with some inference technique. The present implementation, which is based
on the DACE's Matlab toolbox uses the *maximum likelihood estimation* technique.

For a more comprehensive description of the theoretical aspects of Gaussian
Processes for Machine Learning, please refer to the references below:

.. topic:: References:

    * *"DACE, A Matlab Kriging Toolbox"*
      S Lophaven, HB Nielsen, J Sondergaard
      2002,
      <http://www2.imm.dtu.dk/~hbn/dace/>

    * *"Gaussian Processes for Machine Learning"*
      CE Rasmussen, CKI Williams
      MIT Press, 2006 (Ed. T Diettrich)
      <http://www.gaussianprocess.org/gpml/chapters/RW.pdf>

    * *"The design and analysis of computer experiments"*
      TJ Santner, BJ Williams, W Notz
      Springer, 2003
      <http://www.stat.osu.edu/~comp_exp/book.html>

.. correlation_models::

Correlation Models
==================

Common correlation models matches some famous SVM's kernels because they are mostly built on the equivalent
assumptions. They must fulfill Mercer's conditions. Note however, that the choice of the correlation model
should be made in agreement with the known properties of the original experiment from which the observations
come.

* If the original experiment is known to be infinitely differentiable (smooth), then one should use the *squared-exponential correlation model*.
* If it's not, then one should rather use the *exponential correlation model*.
* Note also that there exists a correlation model that takes the degree of derivability as input: this is the Matern correlation model, but it's not implemented here.

For a more detailed discussion on the selection of the appropriate correlation models, dee the book by Rasmussen & Williams in references.

.. regression_models::

Regression Models
=================

Common linear regression models involve zero (constant), first- and second-order polynomials. But one may
specify its own in the form of a Python function that takes the features X as input and that returns a vector
containing the values of the functional set. The only constraint is that the number of functions must not exceed the
number of available observations so that the underlying regression problem is not *under-determined*.


An introductory example
=======================

Say we want to surrogate the function :math:`g(x) = x \sin(x)`. To do so, the function is evaluated onto a
design of experiments. Then, we define a GaussianProcessModel whose regression and correlation
models might be specified using additional kwargs, and ask for the model to be fitted to the data. Depending on the number
of parameters provided at instanciation, the fitting procedure may recourse to maximum likelihood estimation for the parameters
or alternatively it uses the given parameters.



Implementation details
======================

The present implementation is based on a transliteration of the DACE Matlab
toolbox.

.. topic:: References:

    * *"DACE, A Matlab Kriging Toolbox"*
      S Lophaven, HB Nielsen, J Sondergaard
      2002,
      <http://www2.imm.dtu.dk/~hbn/dace/>
