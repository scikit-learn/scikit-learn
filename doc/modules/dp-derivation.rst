.. _dpgmm_derivation:

.. currentmodule:: sklearn.mixture


Variational Gaussian Mixture Models
===================================

The API is identical to that of the :class:`GMM` class, the main
difference being that it offers access to precision matrices as well
as covariance matrices.

The inference algorithm is the one from the following paper:

    * `Variational Inference for Dirichlet Process Mixtures
      <http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.61.4467&rep=rep1&type=pdf>`_
      David Blei, Michael Jordan. Bayesian Analysis, 2006

While this paper presents the parts of the inference algorithm that
are concerned with the structure of the dirichlet process, it does not
go into detail in the mixture modeling part, which can be just as
complex, or even more. For this reason we present here a full
derivation of the inference algorithm and all the update and
lower-bound equations. If you're not interested in learning how to
derive similar algorithms yourself and you're not interested in
changing/debugging the implementation in the scikit this document is
not for you.

The complexity of this implementation is linear in the number of
mixture components and data points. With regards to the
dimensionality, it is linear when using ``spherical`` or ``diag`` and
quadratic/cubic when using ``tied`` or ``full``. For ``spherical`` or ``diag``
it is O(n_states * n_points * dimension) and for ``tied`` or ``full`` it
is O(n_states * n_points * dimension^2 + n_states * dimension^3) (it
is necessary to invert the covariance/precision matrices and compute
its determinant, hence the cubic term).

This implementation is expected to scale at least as well as EM for
the Gaussian mixture.

Update rules for VB inference
==============================

Here the full mathematical derivation of the Variational Bayes update
rules for Gaussian Mixture Models is given. The main parameters of the
model, defined for any class :math:`k \in [1..K]` are the class
proportion :math:`\phi_k`, the mean parameters :math:`\mu_k`, the
covariance parameters :math:`\Sigma_k`, which is characterized by
variational Wishart density, :math:`Wishart(a_k, \mathbf{B_k})`, where
:math:`a` is the degrees of freedom, and :math:`B` is the
scale matrix. Depending on the covariance parametrization,
:math:`B_k` can be a positive scalar, a positive vector or a Symmetric
Positive Definite matrix.


The spherical model
---------------------

The model then is

.. math::

    \begin{array}{rcl}
    \phi_k   &\sim& Beta(1, \alpha_1) \\
    \mu_k   &\sim& Normal(0,  \mathbf{I}) \\
    \sigma_k &\sim& Gamma(1, 1) \\
    z_{i}     &\sim& SBP(\phi) \\
    X_t &\sim& Normal(\mu_{z_i}, \frac{1}{\sigma_{z_i}} \mathbf{I})
    \end{array}

The variational distribution we'll use is

.. math::

    \begin{array}{rcl}
    \phi_k   &\sim& Beta(\gamma_{k,1}, \gamma_{k,2}) \\
    \mu_k   &\sim& Normal(\nu_{\mu_k},  \mathbf{I}) \\
    \sigma_k &\sim& Gamma(a_{k}, b_{k}) \\
    z_{i}     &\sim& Discrete(\nu_{z_i}) \\
    \end{array}


The bound
...........

The variational bound is

.. math::

    \begin{array}{rcl}
    \log P(X) &\ge&
    \sum_k (E_q[\log P(\phi_k)] - E_q[\log Q(\phi_k)]) \\
    &&
    +\sum_k \left( E_q[\log P(\mu_k)] - E_q[\log Q(\mu_k)] \right) \\
    &&
    +\sum_k \left( E_q[\log P(\sigma_k)] - E_q[\log Q(\sigma_k)] \right) \\
    &&
    +\sum_i \left( E_q[\log P(z_i)] - E_q[\log Q(z_i)] \right) \\
    &&
    +\sum_i E_q[\log P(X_t)]
    \end{array}


**The bound for** :math:`\phi_k`

.. math::

    \begin{array}{rcl}
    E_q[\log Beta(1,\alpha)] - E[\log Beta(\gamma_{k,1},\gamma_{k,2})]
    &=&
    \log \Gamma(1+\alpha) - \log \Gamma(\alpha) \\ &&
    +(\alpha-1)(\Psi(\gamma_{k,2})-\Psi(\gamma_{k,1}+\gamma_{k,2})) \\ &&
    - \log \Gamma(\gamma_{k,1}+\gamma_{k,2}) + \log \Gamma(\gamma_{k,1}) +
    \log \Gamma(\gamma_{k,2}) \\ &&
    -
    (\gamma_{k,1}-1)(\Psi(\gamma_{k,1})-\Psi(\gamma_{k,1}+\gamma_{k,2}))
    \\ &&
    -
    (\gamma_{k,2}-1)(\Psi(\gamma_{k,2})-\Psi(\gamma_{k,1}+\gamma_{k,2}))
    \end{array}


**The bound for** :math:`\mu_k`

.. math::

  \begin{array}{rcl}
  && E_q[\log P(\mu_k)] - E_q[\log Q(\mu_k)] \\
  &=&
  \int\!d\mu_f q(\mu_f) \log P(\mu_f)
  - \int\!d\mu_f q(\mu_f) \log Q(\mu_f)  \\
  &=&
  - \frac{D}{2}\log 2\pi - \frac{1}{2} ||\nu_{\mu_k}||^2 - \frac{D}{2}
  + \frac{D}{2} \log 2\pi e
  \end{array}


**The bound for** :math:`\sigma_k`

Here I'll use the inverse scale parametrization of the gamma
distribution.

.. math::

  \begin{array}{rcl}
  && E_q[\log P(\sigma_k)] - E_q [\log Q(\sigma_k)] \\ &=&
  \log \Gamma (a_k) - (a_k-1)\Psi(a_k) -\log b_k + a_k - \frac{a_k}{b_k}
  \end{array}


**The bound for z**

.. math::

  \begin{array}{rcl}
  && E_q[\log P(z)] - E_q[\log Q(z)] \\
  &=&
  \sum_{k} \left(
       \left(\sum_{j=k+1}^K  \nu_{z_{i,j}}\right)(\Psi(\gamma_{k,2})-\Psi(\gamma_{k,1}+\gamma_{k,2}))
   +  \nu_{z_{i,k}}(\Psi(\gamma_{k,1})-\Psi(\gamma_{k,1}+\gamma_{k,2}))
   - \log \nu_{z_{i,k}} \right)
  \end{array}


**The bound for** :math:`X`

Recall that there is no need for a :math:`Q(X)` so this bound is just

.. math::

    \begin{array}{rcl}
    E_q[\log P(X_i)] &=& \sum_k \nu_{z_k} \left( - \frac{D}{2}\log 2\pi
    +\frac{D}{2} (\Psi(a_k) - \log(b_k))
    -\frac{a_k}{2b_k} (||X_i - \nu_{\mu_k}||^2+D) - \log 2 \pi e  \right)
    \end{array}


For simplicity I'll later call the term inside the parenthesis :math:`E_q[\log P(X_i|z_i=k)]`

The updates
............

**Updating** :math:`\gamma`

.. math::

  \begin{array}{rcl}
  \gamma_{k,1} &=& 1+\sum_i \nu_{z_{i,k}} \\
  \gamma_{k,2} &=& \alpha + \sum_i \sum_{j > k} \nu_{z_{i,j}}.
  \end{array}


**Updating** :math:`\mu`

The updates for mu essentially are just weighted expectations of
:math:`X` regularized by the prior. We can see this by taking the
gradient of the bound with regards to :math:`\nu_{\mu}` and setting it to zero.
The gradient is

.. math::

  \nabla L = -\nu_{\mu_k} + \sum_i \frac{\nu_{z_{i,k}}b_k}{a_k}(X_i + -\nu_{\mu})


so the update is

.. math::

    \nu_{\mu_k} = \frac{\sum_i \frac{\nu_{z_{i,k}}b_k}{a_k}X_i}{1+\sum_i \frac{\nu_{z_{i,k}}b_k}{a_k}}



**Updating** :math:`a` **and** :math:`b`


For some odd reason it doesn't really work when you derive the updates
for a and b using the gradients of the lower bound (terms involving the
:math:`\Psi'` function show up and :math:`a` is hard to isolate).
However, we can use the other formula,

.. math::

  \log Q(\sigma_k) = E_{v \ne \sigma_k}[\log P] + const


All the terms not involving :math:`\sigma_k` get folded over into the
constant and we get two terms: the prior and the probability of
:math:`X`. This gives us

.. math::

   \log Q(\sigma_k) = -\sigma_k  + \frac{D}{2} \sum_i \nu_{z_{i,k}}\log \sigma_k  - \frac{\sigma_k}{2}\sum_i \nu_{z_{i,k}} (||X_i-\mu_k||^2 + D)


This is the log of a gamma distribution, with :math:`a_k = 1+\frac{D}{2}\sum_i \nu_{z_{i,k}}` and

.. math::

  b_k = 1 + \frac{1}{2}\sum_i \nu_{z_{i,k}} (||X_i-\mu_k||^2 + D).


You can verify this by normalizing the previous term.

**Updating** :math:`z`

.. math::

   \log \nu_{z_{i,k}} \propto \Psi(\gamma_{k,1}) -
   \Psi(\gamma_{k,1} + \gamma_{k,2}) + E_q[\log P(X_i|z_i=k)] +
   \sum_{j < k} \left (\Psi(\gamma_{j,2}) -
   \Psi(\gamma_{j,1}+\gamma_{j,2})\right).


The diagonal model
--------------------


The model then is

.. math::

  \begin{array}{rcl}
  \phi_k   &\sim& Beta(1, \alpha_1) \\
  \mu_k   &\sim& Normal(0,  \mathbf{I}) \\
  \sigma_{k,d} &\sim& Gamma(1, 1) \\
  z_{i}     &\sim& SBP(\phi) \\
  X_t &\sim& Normal(\mu_{z_i}, \bm{\sigma_{z_i}}^{-1})
  \end{array}

Tha variational distribution we'll use is

.. math::

  \begin{array}{rcl}
  \phi_k   &\sim& Beta(\gamma_{k,1}, \gamma_{k,2}) \\
  \mu_k   &\sim& Normal(\nu_{\mu_k},  \mathbf{I}) \\
  \sigma_{k,d} &\sim& Gamma(a_{k,d}, b_{k,d}) \\
  z_{i}     &\sim& Discrete(\nu_{z_i}) \\
  \end{array}

The lower bound
...................

The changes in this lower bound from the previous model are in the
distributions of :math:`\sigma` (as there are a lot more :math:`\sigma` s
now) and :math:`X`.

The bound for :math:`\sigma_{k,d}` is the same bound for :math:`\sigma_k` and can
be safely omitted.

**The bound for** :math:`X` :

The main difference here is that the precision matrix :math:`\bm{\sigma_k}`
scales the norm, so we have an extra term after computing the
expectation of :math:`\mu_k^T\bm{\sigma_k}\mu_k`, which is
:math:`\nu_{\mu_k}^T\bm{\sigma_k}\nu_{\mu_k} + \sum_d \sigma_{k,d}`. We then
have

.. math::

  \begin{array}{rcl}
  E_q[\log P(X_i)] &=& \sum_k \nu_{z_k} \Big( - \frac{D}{2}\log 2\pi
  +\frac{1}{2}\sum_d (\Psi(a_{k,d}) - \log(b_{k,d})) \\
  &&
  -\frac{1}{2}((X_i - \nu_{\mu_k})^T\bm{\frac{a_k}{b_k}}(X_i - \nu_{\mu_k})+ \sum_d \sigma_{k,d})- \log 2 \pi e  \Big)
  \end{array}


The updates
............

The updates only chance for :math:`\mu` (to weight them with the new
:math:`\sigma`), :math:`z` (but the change is all folded into the
:math:`E_q[P(X_i|z_i=k)]` term), and the :math:`a` and :math:`b` variables themselves.

**The update for** :math:`\mu`

.. math::

   \nu_{\mu_k} = \left(\mathbf{I}+\sum_i \frac{\nu_{z_{i,k}}\mathbf{b_k}}{\mathbf{a_k}}\right)^{-1}\left(\sum_i \frac{\nu_{z_{i,k}}b_k}{a_k}X_i\right)


**The updates for a and b**

Here we'll do something very similar to the spheric model. The main
difference is that now each :math:`\sigma_{k,d}` controls only one dimension
of the bound:

.. math::

  \log Q(\sigma_{k,d}) = -\sigma_{k,d} + \sum_i \nu_{z_{i,k}}\frac{1}{2}\log \sigma_{k,d}
  - \frac{\sigma_{k,d}}{2}\sum_i \nu_{z_{i,k}} ((X_{i,d}-\mu_{k,d})^2 + 1)


Hence

.. math::

  a_{k,d} = 1 + \frac{1}{2} \sum_i \nu_{z_{i,k}}

.. math::

  b_{k,d} = 1 + \frac{1}{2} \sum_i \nu_{z_{i,k}}((X_{i,d}-\mu_{k,d})^2 + 1)


The tied model
----------------

The model then is

.. math::

  \begin{array}{rcl}
  \phi_k   &\sim& Beta(1, \alpha_1) \\
  \mu_k   &\sim& Normal(0,  \mathbf{I}) \\
  \Sigma &\sim& Wishart(D, \mathbf{I}) \\
  z_{i}     &\sim& SBP(\phi) \\
  X_t &\sim& Normal(\mu_{z_i},  \Sigma^{-1})
  \end{array}

Tha variational distribution we'll use is

.. math::

  \begin{array}{rcl}
  \phi_k   &\sim& Beta(\gamma_{k,1}, \gamma_{k,2}) \\
  \mu_k   &\sim& Normal(\nu_{\mu_k},  \mathbf{I}) \\
  \Sigma &\sim& Wishart(a, \mathbf{B}) \\
  z_{i}     &\sim& Discrete(\nu_{z_i}) \\
  \end{array}

The lower bound
..................

There are two changes in the lower-bound: for :math:`\Sigma` and for :math:`X`.

**The bound for** :math:`\Sigma`

.. math::

  \begin{array}{rcl}
  \frac{D^2}{2}\log 2  + \sum_d \log \Gamma(\frac{D+1-d}{2}) \\
  - \frac{aD}{2}\log 2 + \frac{a}{2} \log |\mathbf{B}| + \sum_d \log \Gamma(\frac{a+1-d}{2}) \\
  + \frac{a-D}{2}\left(\sum_d \Psi\left(\frac{a+1-d}{2}\right)
  + D \log 2 + \log |\mathbf{B}|\right) \\
  + \frac{1}{2} a \mathbf{tr}[\mathbf{B}-\mathbf{I}]
  \end{array}


**The bound for X**

.. math::

   \begin{array}{rcl}
   E_q[\log P(X_i)] &=& \sum_k \nu_{z_k} \Big( - \frac{D}{2}\log 2\pi
   +\frac{1}{2}\left(\sum_d \Psi\left(\frac{a+1-d}{2}\right)
   + D \log 2 + \log |\mathbf{B}|\right) \\
   &&
   -\frac{1}{2}((X_i - \nu_{\mu_k})a\mathbf{B}(X_i - \nu_{\mu_k})+ a\mathbf{tr}(\mathbf{B}))- \log 2 \pi e  \Big)
   \end{array}

The updates
.............

As in the last setting, what changes are the trivial update for :math:`z`,
the update for :math:`\mu` and the update for :math:`a` and :math:`\mathbf{B}`.

**The update for** :math:`\mu`

.. math::

    \nu_{\mu_k} = \left(\mathbf{I}+ a\mathbf{B}\sum_i \nu_{z_{i,k}}\right)^{-1}
    \left(a\mathbf{B}\sum_i \nu_{z_{i,k}} X_i\right)

**The update for** :math:`a` **and** :math:`B`

As this distribution is far too complicated I'm not even going to try
going at it the gradient way.

.. math::

   \log Q(\Sigma) = +\frac{1}{2}\log |\Sigma| - \frac{1}{2} \mathbf{tr}[\Sigma]
   + \sum_i \sum_k \nu_{z_{i,k}} \left( +\frac{1}{2}\log |\Sigma| - \frac{1}{2}((X_i-\nu_{\mu_k})^T\Sigma(X_i-\nu_{\mu_k})+\mathbf{tr}[\Sigma]) \right)

which non-trivially (seeing that the quadratic form with :math:`\Sigma` in
the middle can be expressed as the trace of something) reduces to

.. math::

   \log Q(\Sigma) = +\frac{1}{2}\log |\Sigma| - \frac{1}{2} \mathbf{tr}[\Sigma]
   + \sum_i \sum_k \nu_{z_{i,k}} \left( +\frac{1}{2}\log |\Sigma| - \frac{1}{2}(\mathbf{tr}[(X_i-\nu_{\mu_k})(X_i-\nu_{\mu_k})^T\Sigma]+\mathbf{tr}[I \Sigma]) \right)

hence this (with a bit of squinting) looks like a wishart with parameters

.. math::

   a = 2 + D + T


and

.. math::

   \mathbf{B} = \left(\mathbf{I} + \sum_i \sum_k \nu_{z_{i,k}}(X_i-\nu_{\mu_k})(X_i-\nu_{\mu_k})^T\right)^{-1}


The full model
----------------

 The model then is

.. math::

  \begin{array}{rcl}
  \phi_k   &\sim& Beta(1, \alpha_1) \\
  \mu_k   &\sim& Normal(0,  \mathbf{I}) \\
  \Sigma_k &\sim& Wishart(D, \mathbf{I}) \\
  z_{i}     &\sim& SBP(\phi) \\
  X_t &\sim& Normal(\mu_{z_i},  \Sigma_{z,i}^{-1})
  \end{array}


The variational distribution we'll use is

.. math::

  \begin{array}{rcl}
  \phi_k   &\sim& Beta(\gamma_{k,1}, \gamma_{k,2}) \\
  \mu_k   &\sim& Normal(\nu_{\mu_k},  \mathbf{I}) \\
  \Sigma_k &\sim& Wishart(a_k, \mathbf{B_k}) \\
  z_{i}     &\sim& Discrete(\nu_{z_i}) \\
  \end{array}

The lower bound
.................

All that changes in this lower bound in comparison to the previous one
is that there are K priors on different :math:`\Sigma` precision matrices
and there are the correct indices on the bound for X.

The updates
..............

All that changes in the updates is that the update for mu uses only
the proper sigma and the updates for a and B don't have a sum over K, so

.. math::

         \nu_{\mu_k} = \left(\mathbf{I}+ a_k\mathbf{B_k}\sum_i \nu_{z_{i,k}}\right)^{-1}
         \left(a_k\mathbf{B_k}\sum_i \nu_{z_{i,k}} X_i\right)

.. math::

    a_k = 2 + D + \sum_i \nu_{z_{i,k}}

and

.. math::

       \mathbf{B} = \left(\left(\sum_i\nu_{z_{i,k}}+1\right)\mathbf{I} + \sum_i  \nu_{z_{i,k}}(X_i-\nu_{\mu_k})(X_i-\nu_{\mu_k})^T\right)^{-1}

