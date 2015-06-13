.. _mixture_derivation:

.. currentmodule:: sklearn.mixture


Bayesian Gaussian mixture with variational inference
====================================================
Here we briefly represent the updating equations of Bayesian Gaussian mixture
with variational inference. The inference algorithm is the one from the book
Pattern Recognition and Machine Learning:

    * `Pattern Recognition and Machine Learning
      <http://research.microsoft.com/en-us/um/people/cmbishop/prml/>`_
      Christopher M. Bishop, 2007


We denote latent variable :math:`\mathbf{z}_n` (1-of-:math:`K` encoding),
the observed data set :math:`\mathbf{X}=\{\mathbf{x}_1, \ldots, \mathbf{x}_N \}`,
:math:`\mathbf{Z} = \{\mathbf{z}_1, \ldots, \mathbf{z}_N \}`.
The joint probability of Bayesian mixture model is

.. math::

    p(\mathbf{X}, \mathbf{Z}, \boldsymbol{\pi}, \boldsymbol{\mu}, \boldsymbol{\Lambda}) =
    p(\mathbf{X} | \mathbf{Z}, \boldsymbol{\mu}, \boldsymbol{\Lambda})
    p(\mathbf{Z} | \boldsymbol{\pi})
    p(\boldsymbol{\pi})
    p(\boldsymbol{\mu}|\boldsymbol{\Lambda})
    p(\boldsymbol{\Lambda})

and the following variables and distributions.

.. math::

    \begin{array}{rl}
        N_k & = \sum_{n=1}^N r_{nk} \\
        \bar{\mathbf{x}}_k & = \frac{1}{N_k} \sum_{n=1}^Nr_{nk}\mathbf{x}_n \\
        \mathbf{S}_{k} & = \frac{1}{N_k} \sum_{n=1}^N r_{nk} (\mathbf{x}_{n} - \bar{\mathbf{x}}_k)(\mathbf{x}_{n} - \bar{\mathbf{x}}_k)^\top
    \end{array}

.. math::

    \begin{array}{rl}
        \mathcal{W}(\boldsymbol{\Lambda}|\mathbf{W}, \nu) & = B(\mathbf{W}, \nu)|\boldsymbol{\Lambda}|^{(\nu-D-1)/2}
          \exp\left( -\frac{1}{2} \textrm{Tr}(\mathbf{W}^{-1}\boldsymbol{\Lambda})\right) \\
        \textrm{Gam} (\tau|a, b) & = G(a, b) \tau^{a-1} e^{- \frac{\tau}{b}} \\
        \textrm{Dir} (\boldsymbol{\mu}|\boldsymbol{\alpha}) &= C(\boldsymbol{\alpha}) \prod_{k=1}^K \mu_k^{\alpha_k - 1} \\
        \textrm{Beta}(\mu|a, b) & = F(a, b)\mu^{a-1}(1-\mu)^{b-1}
    \end{array}

The terms in the distributions :math:`\textrm{B}`, :math:`\textrm{G}`, :math:`\textrm{C}`, :math:`\textrm{F}`
are normalization terms.

Full Precision
--------------
When the precision is ``full`` and ``tied``, we use Wishart distribution for the
precision. So in the joint probability, each term has the following form:

.. math::
    \begin{array}{rl}
        p(\mathbf{X} | \mathbf{Z}, \boldsymbol{\mu}, \boldsymbol{\Lambda})
            & =  \prod_{n=1}^N\prod_{k=1}^K \mathcal{N} (\mathbf{x}_n | \boldsymbol{\mu}_k, \boldsymbol{\Lambda}_k^{-1})^{z_{nk}} \\
        p(\mathbf{Z} | \boldsymbol{\pi})
            & =  \prod_{n=1}^N \prod_{k=1}^K \pi_k^{z_{nk}} \\
        p(\boldsymbol{\pi})
            & =  \textrm{Dir}(\boldsymbol{\pi}| \boldsymbol{\alpha}_0) = C(\boldsymbol{\alpha}_0) \prod_{k=1}^K \pi_k^{\alpha_0-1} \\
        p(\boldsymbol{\mu}, \boldsymbol{\Lambda})
            & =  p(\boldsymbol{\mu}|\boldsymbol{\Lambda})p(\boldsymbol{\Lambda})  \\
            & = \prod_{k=1}^K \mathcal{N}(\boldsymbol{\mu}_k|\mathbf{m}_0, (\beta_0\boldsymbol{\Lambda}_k)^{-1})
                          \mathcal{W}(\boldsymbol{\Lambda}_k| \mathbf{W}_0, \nu_0) \\
        q(\mathbf{Z}, \boldsymbol{\pi}, \boldsymbol{\mu}, \boldsymbol{\Lambda}) & =
        q(\mathbf{Z})q(\boldsymbol{\pi}, \boldsymbol{\mu}, \boldsymbol{\Lambda}) \\
        q(\boldsymbol{\pi}, \boldsymbol{\mu}, \boldsymbol{\Lambda}) & =
        q(\boldsymbol{\pi})\prod_{k=1}^K q(\boldsymbol{\mu}_k, \boldsymbol{\Lambda}_k)
    \end{array}

E step
......

The responsibility of each data point :math:`r_{nk}`  is the normalized :math:`\rho_{nk}`.

.. math::

    \begin{array}{rl}
        q^*(\mathbf{Z}) & = \prod_{n=1}^N\prod_{k=1}^K r_{nk}^{z_{nk}} \\
        r_{nk} & = \frac{\rho_{nk}}{\sum_{i=1}^K \rho_{ni}}
    \end{array}

:math:`\rho_{nk}` is computed from the following equation.

.. math::
    \ln \rho_{nk} = \mathbb{E}[\ln \pi_k] +
                    \frac{1}{2} \mathbb{E}[\ln |\boldsymbol{\Lambda}_k|] -
                    \frac{D}{2} \ln(2\pi) -
                    \frac{1}{2} \mathbb{E}_{\boldsymbol{\mu}_k, \boldsymbol{\Lambda}_k}
                                \left[
                                  (\mathbf{x}_n - \boldsymbol{\mu}_k)^\top \boldsymbol{\Lambda}_k (\mathbf{x}_n - \boldsymbol{\mu}_k)
                                \right]

.. math::

    \begin{array}{rl}
        \ln \widetilde{\pi}_k = \mathbb{E}[\ln \pi_k] & = \psi(\alpha_k) - \psi(\hat{\alpha}) \\
        \ln \widetilde{\Lambda}_k = \mathbb{E}[\ln |\boldsymbol{\Lambda}_k|]
            & = \sum_{i=1}^D \psi \left( \frac{\nu_k + 1 - i}{2} \right) +  D \ln 2 + \ln |\mathbf{W}_k|
    \end{array}


.. math::

    \begin{array}{rl}
    \mathbb{E}_{\boldsymbol{\mu}_k, \boldsymbol{\Lambda}_k}
        \left[
            (\mathbf{x}_n - \boldsymbol{\mu}_k)^\top \boldsymbol{\Lambda}_k (\mathbf{x}_n - \boldsymbol{\mu}_k)
        \right]
     = D \beta_k^{-1} + \nu_k (\mathbf{x}_n - \mathbf{m}_k)^\top \mathbf{W}_k (\mathbf{x}_n - \mathbf{m}_k)
    \end{array}

M step
......

In M-step, the parameters of :math:`\boldsymbol{\pi}`, :math:`\boldsymbol{\mu}_k`
and :math:`\boldsymbol{\Lambda}_k` distributions should be updated.

.. math::

    \begin{array}{rl}
        q^*(\boldsymbol{\pi}) & = \textrm{Dir} (\boldsymbol{\pi}|\boldsymbol{\alpha}) \\
        \alpha_k & = \alpha_0 + N_k
    \end{array}


.. math::

    \begin{array}{rl}
        q^*(\boldsymbol{\mu}_k| \boldsymbol{\Lambda}_k) & =
            \mathcal{N} (\boldsymbol{\mu}_k | \mathbf{m}_k, (\beta_k \boldsymbol{\Lambda}_k)^{-1}) \\
        \beta_k & = \beta_0 + N_k \\
        \mathbf{m}_k & = \frac{1}{\beta_k}(\beta_0 \mathbf{m}_0 + N_k \bar{\mathbf{x}}_k)
    \end{array}


.. math::

    \begin{array}{rl}
        q^*(\boldsymbol{\Lambda}_k) &= \mathcal{W}(\boldsymbol{\Lambda}_k | \mathbf{W}_k,  \nu_k) \\
        \nu_k & = \nu_0 + N_k \\
        \mathbf{W}_k^{-1} & = \mathbf{W}_0^{-1} + N_k \mathbf{S}_k +
        \frac{N_k \beta_0}{\beta_k}(\bar{\mathbf{x}}_k - \mathbf{m}_0)(\bar{\mathbf{x}}_k - \mathbf{m}_0)^\top
    \end{array}


Tied Precision
--------------

The updating functions of Bayesian mixture model with tied precision are the
same as those of full precision, except that

.. math::

    \begin{array}{rl}
        \nu & = \nu_0 + \frac{N}{K} \\
        \mathbf{C}_k & = \mathbf{W}_0^{-1} + N_k \mathbf{S}_k +
        \frac{N_k \beta_0}{\beta_k}(\bar{\mathbf{x}}_k - \mathbf{m}_0)(\bar{\mathbf{x}}_k - \mathbf{m}_0)^\top \\
        \mathbf{W}^{-1} & = \frac{1}{K}\sum_{k=1}^K \mathbf{C}_k
    \end{array}


Diag Precision
--------------
When precision matrices are diagonal for all components, we assume that in each dimension the diagonal variables
of precision matrices follow Gamma distribution,
:math:`\boldsymbol{\tau} = \{\boldsymbol{\tau}_1, \boldsymbol{\tau}_2, \ldots, \boldsymbol{\tau}_K \}`,
:math:`\boldsymbol{\tau}_k = \textrm{diag} \{ \tau_{k1}, \tau_{k2}, \ldots, \tau_{kD}\}`.
It should be noticed that for each component :math:`k`, there are :math:`D` Gaussian-Gamma distribution.
The :math:`K \times D` distributions have common prior parameters :math:`a_0` and :math:`\mathbf{b}_0`.
Therefore the :math:`D` posterior distributions for each component have
the posterior parameters :math:`a_k` and :math:`\mathbf{b}_k`. :math:`\mathbf{b}_0 = \textrm{diag} \{b_{01}, b_{02}, \ldots, b_{0D} \}`.
:math:`\mathbf{b}_k = \textrm{diag} \{b_{k1}, b_{k2}, \ldots, b_{kD} \}`.
In other words, these :math:`D` distributions have their own :math:`b_{kd}`, but share :math:`a_k`.

.. math::

    \begin{array}{rl}
      p(\mathbf{X} | \mathbf{Z}, \boldsymbol{\mu}, \boldsymbol{\Lambda})
        & = \prod_{n=1}^N \prod_{k=1}^K \prod_{d=1}^D
            \mathcal{N} (x_{nd} | \mu_{kd}, \tau_{kd}^{-1})^{z_{nk}} \\
      p(\mathbf{Z} | \boldsymbol{\pi})
        & =  \prod_{n=1}^N \prod_{k=1}^K \pi_k^{z_{nk}} \\
      p(\boldsymbol{\pi})
        & =  \textrm{Dir}(\boldsymbol{\pi}| \boldsymbol{\alpha}_0) \prod_{k=1}^K \pi_k^{\alpha_0-1} \\
      p(\boldsymbol{\mu}, \boldsymbol{\tau})
        & =  p(\boldsymbol{\mu}|\boldsymbol{\tau})p(\boldsymbol{\tau})  \\
        & = \prod_{k=1}^K \prod_{d=1}^D
              \mathcal{N}(\mu_{kd}|m_{0d}, (\lambda_0 \tau_{kd})^{-1})
                          \textrm{Gam} (\tau_{kd}| a_0, b_{0d}) \\
      q(\mathbf{Z}, \boldsymbol{\pi}, \boldsymbol{\mu}, \boldsymbol{\tau}) & =
      q(\mathbf{Z})q(\boldsymbol{\pi}, \boldsymbol{\mu}, \boldsymbol{\tau}) \\
      q(\boldsymbol{\pi}, \boldsymbol{\mu}, \boldsymbol{\tau}) & =
      q(\boldsymbol{\pi})\prod_{k=1}^K q(\boldsymbol{\mu}_k, \boldsymbol{\tau}_k)
    \end{array}


The updating functions of Bayesian mixture model with ``diag`` precision are the
same as those of full precision, except that in M-step,

.. math::

    \begin{array}{rl}
        a_k &= a_0 + \frac{N_k}{2} \\
        q^*(\boldsymbol{\tau}_k) &= \prod_{d=1}^D \textrm{Gam} (\tau_{kd} | a_k, b_{kd}) \\
        a_k & = a_0 + \frac{N_k}{2} \\
        b_{kd}^{-1} &= b_{0d}^{-1} + \frac{1}{2} \textrm{diag}
            \left[
                N_k \mathbf{S}_k +
                \frac{N_k \lambda_0}{\lambda_k}(\bar{\mathbf{x}}_k - \mathbf{m}_0)(\bar{\mathbf{x}}_k - \mathbf{m}_0)^\top
            \right]_{d}
    \end{array}

Spherical Precision
-------------------

The M-step for Bayesian Gaussian mixture with ``spherical`` precision is

.. math::

    \begin{array}{rl}
        q^*(\boldsymbol{\tau}_k) &= \prod_{d=1}^D \textrm{Gam} (\tau_{kd} | a_k, b_{kd}) \\
        a_k & = a_0 + \frac{N_k}{2} \\
        b_{k}^{-1} &= b_{0}^{-1} + \frac{1}{2D} \textrm{Tr}
            \left[
                N_k \mathbf{S}_k +
                \frac{N_k \lambda_0}{\lambda_k}(\bar{\mathbf{x}}_k - \mathbf{m}_0)(\bar{\mathbf{x}}_k - \mathbf{m}_0)^\top
            \right]
    \end{array}


Lower bound (ELBO)
------------------
The lower bound of the probability distribution is given by

.. math::

    \begin{array}{rl}
        \mathcal{L}
      & =  \sum_{\mathbf{Z}} \iiint q(\mathbf{Z}, \boldsymbol{\pi}, \boldsymbol{\mu}, \boldsymbol{\Lambda})
            \ln \left\{
            \frac{p(\mathbf{X}, \mathbf{Z}, \boldsymbol{\pi}, \boldsymbol{\mu}, \boldsymbol{\Lambda})}
                 {q(\mathbf{Z}, \boldsymbol{\pi}, \boldsymbol{\mu}, \boldsymbol{\Lambda})}\right\}
            \,\mathrm{d} \boldsymbol{\pi} \,\mathrm{d} \boldsymbol{\mu} \,\mathrm{d} \boldsymbol{\Lambda} \\
      & =  \mathbb{E}[\ln p(\mathbf{X}|\mathbf{Z},\boldsymbol{\mu}, \boldsymbol{\Lambda})]
            + \mathbb{E} [\ln p(\mathbf{Z}|\boldsymbol{\pi})]
            + \mathbb{E} [\ln p(\boldsymbol{\pi})]
            + \mathbb{E} [\ln p(\boldsymbol{\mu}, \boldsymbol{\Lambda})] \\
       &  - \mathbb{E} [\ln q(\mathbf{Z})]
            - \mathbb{E} [\ln q(\boldsymbol{\pi})]
            - \mathbb{E} [\ln q(\boldsymbol{\mu}, \boldsymbol{\Lambda})]
    \end{array}

We need to compute 7 terms, but the first two terms can be computed from the responsibility, the product
of the E-step.

.. math::

    \begin{array}{rl}
        \mathbb{E}[\ln p(\mathbf{X}|\mathbf{Z},\boldsymbol{\mu}, \boldsymbol{\Lambda})]
            & =  \frac{1}{2} \mathbb{E}
            \left[
                \sum_{n=1}^N \sum_{k=1}^K z_{nk}
                \left\{
                  \ln | \boldsymbol{\Lambda}_k| - D \ln (2\pi)
                  - (\mathbf{x}_n - \boldsymbol{\mu}_k)^\top \boldsymbol{\Lambda}_k (\mathbf{x}_n - \boldsymbol{\mu}_k)
                \right\}
           \right] \\
      & = \frac{1}{2} \sum_{k=1}^K N_k \left\{ \ln \widetilde{\Lambda}_k - D \ln(2\pi) \right\}
            - \sum_{n=1}^N \sum_{k=1}^K r_{nk} [ D\beta_k^{-1} + \nu_k (\mathbf{x}_n - \mathbf{m}_k)^\top \mathbf{W}_k
            (\mathbf{x}_n - \mathbf{m}_k)] \\
      \label{eq:lb_px}
      & = \frac{1}{2} \sum_{k=1}^K N_k \left\{ \ln \widetilde{\Lambda}_k - D \ln(2\pi)
            - D\beta_k^{-1}
            - \nu_k \textrm{Tr}[\mathbf{S}_k \mathbf{W}_k]
            - \nu_k (\bar{\mathbf{x}}_k - \mathbf{m}_k)^\top \mathbf{W}_k (\bar{\mathbf{x}}_n - \mathbf{m}_k)]
            \right\}
    \end{array}

.. math::

    \begin{array}{rl}
        \mathbb{E}[\ln p(\mathbf{Z}|\boldsymbol{\pi})]
        & = \sum_{n=1}^N \sum_{k=1}^K r_{nk} \ln \widetilde{\pi}_k \\
        \mathbb{E}[\ln p(\boldsymbol{\pi})]
        & = \ln C(\boldsymbol{\alpha}_0) + (\alpha_0 - 1) \sum_{k=1}^K \ln \widetilde{\pi}_k
    \end{array}

.. math::

    \begin{array}{rl}
        \mathbb{E} [\ln p(\boldsymbol{\mu}, \boldsymbol{\Lambda})]
        & = \sum_{k=1}^K \mathbb{E}
            \left[ \ln \mathcal{N}(\boldsymbol{\mu}_k | \mathbf{m}_0, (\beta_0 \boldsymbol{\Lambda}_k)^{-1}) +
                 \ln \mathcal{W}(\boldsymbol{\Lambda}_k | \mathbf{W}_0, \nu_0 )
          \right] \\
        & = \frac{1}{2} \sum_{k=1}^K
          \left\{
              D \ln \frac{\beta_0}{2\pi}
              + \ln \widetilde{\Lambda}_k
              - \frac{D \beta_0}{\beta_k}
              - \beta_0 \nu_k (\mathbf{m}_k - \mathbf{m}_0)^\top \mathbf{W}_k (\mathbf{m}_k - \mathbf{m}_0)
          \right\} \\
        &   + K \ln B(\mathbf{W}_0, \nu_0) + \frac{\nu_0-D-1}{2} \sum_{k=1}^K \ln \widetilde{\Lambda}_k
              - \frac{1}{2} \sum_{k=1}^K \nu_k \textrm{Tr} (\mathbf{W}_0^{-1}\mathbf{W}_k)
    \end{array}

.. math::

    \begin{array}{rl}
        \mathbb{E}[\ln q(\mathbf{Z})] & = \sum_{n=1}^N \sum_{k=1}^K r_{nk} \ln r_{nk} \\
        \mathbb{E}[\ln q(\boldsymbol{\pi})] & = \sum_{k=1}^K(\alpha_k - 1) \ln \tilde{\pi}_k + \ln C(\boldsymbol{\alpha})
    \end{array}

.. math::

    \begin{array}{rl}
    \mathbb{E}[\ln q(\boldsymbol{\mu}, \boldsymbol{\Lambda})]
      & = \sum_{k=1}^K \mathbb{E}[
          \ln \mathcal{N}(\boldsymbol{\mu}_k | \mathbf{m}_k, (\beta_k \boldsymbol{\Lambda}_k)^{-1}) +
          \ln \mathcal{W}(\boldsymbol{\Lambda}_k | \mathbf{W}_k, \nu_k) ] \\
      & = \sum_{k=1}^K
          \left\{
            \frac{1}{2} \ln \widetilde{\Lambda}_k
            + \frac{D}{2} \ln \frac{\beta_k}{2 \pi}
            - \frac{1}{2} \mathbb{E}[(\boldsymbol{\mu}_k - \mathbf{m}_k)^\top \beta_k\boldsymbol{\Lambda}_k (\boldsymbol{\mu}_k - \mathbf{m}_k)]
            - \mathrm{H} [q(\boldsymbol{\Lambda}_k)]
          \right\} \\
      & = \sum_{k=1}^K
          \left\{
            \frac{1}{2} \ln \widetilde{\Lambda}_k
            + \frac{D}{2} \ln \frac{\beta_k}{2 \pi}
            - \frac{D}{2}
            - \mathrm{H} [ q( \boldsymbol{\Lambda}_k )]
          \right\}
    \end{array}


The entropy of Wishart distribution is defined as

.. math::

  \mathrm{H}[\boldsymbol{\Lambda}] =
      -\ln B(\mathbf{W}, \nu)
      - \frac{(\nu - D - 1)}{2} \mathbb{E}[\ln | \boldsymbol{\Lambda}|]
      + \frac{\nu D}{2}


Dirichlet Process Gaussian Mixture
----------------------------------

Dirichlet process Gaussian mixture model is an extension of Bayesian GMM.

    * `Variational Inference for Dirichlet Process Mixtures
      <http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.61.4467&rep=rep1&type=pdf>`_
      David Blei, Michael Jordan. Bayesian Analysis, 2006

The author used the truncated stick-breaking approach on the
variational distribution to approximate the model.
The stick-breaking approach is like randomly breaking a unit-length stick infinite times.
Suppose we have two infinite
collections of independent random variables :math:`\mathbf{V}` and :math:`\boldsymbol{\eta}^*`,
the stick-breaking representation of Dirichlet process is

.. math::

    \begin{array}{rl}
        V_i & \sim \textrm{Beta}(V_i | 1, \gamma_0) \\
        \eta_i^* & \sim G_0 \\
        \pi_i & = V_i \prod_{j=1}^{i-1}(1-V_j)  \\
        G &= \sum_{i=1}^{\infty} \pi_i \delta_{\eta_i^*}
    \end{array}

Here :math:`G_0` could be the Gaussian distribution, and
:math:`G` is the Dirichlet process Gaussian mixture.
Then the variational inference is applied on the stick-breaking representation of the DPGM.
It should be noticed the coordinate ascent algorithm for mean-field variational inference in the Blei's paper
is exactly the same as the variational inference in PRML.
The stick-breaking representation of the variational distributions is truncated.

The joint probability is

.. math::
    p(\mathbf{X}, \mathbf{Z}, \boldsymbol{\pi}, \boldsymbol{\mu}, \boldsymbol{\Lambda}) =
    p(\mathbf{X} | \mathbf{Z}, \boldsymbol{\mu}, \boldsymbol{\Lambda})
    p(\mathbf{Z} | \mathbf{V})
    p(\mathbf{V})
    p(\boldsymbol{\mu}|\boldsymbol{\Lambda})
    p(\boldsymbol{\Lambda})

The mixture coefficient :math:`\boldsymbol{\pi}` is replaced by an infinity collection :math:`\mathbf{V}`.
The terms are

.. math::
    \begin{array}{rl}
    p(\mathbf{X} | \mathbf{Z}, \boldsymbol{\mu}, \boldsymbol{\Lambda}) & = \prod_{n=1}^N\prod_{k=1}^\infty \mathcal{N} (\mathbf{x}_n | \boldsymbol{\mu}_k, \boldsymbol{\Lambda}_k^{-1})^{z_{nk}} \\
    p(\mathbf{Z} | \mathbf{V}) & = \prod_{n=1}^N \prod_{k=1}^\infty \left[ V_k \prod_{i=1}^{k-1}(1 - V_i)\right]^{z_{nk}} \\
    p(\mathbf{V}) & = \prod_{k=1}^\infty \textrm{Beta}(V_i|1, \gamma_0) \\
    p(\boldsymbol{\mu}, \boldsymbol{\Lambda}) & = p(\boldsymbol{\mu}|\boldsymbol{\Lambda})p(\boldsymbol{\Lambda})  \\
    & = \prod_{k=1}^\infty \mathcal{N}(\boldsymbol{\mu}_k| \mathbf{m}_0, (\beta_0\boldsymbol{\Lambda}_k)^{-1}) \mathcal{W}(\boldsymbol{\Lambda}_k| \mathbf{W}_0, \nu_0)
    \end{array}

The truncated variational distributions are, (we use subscript :math:`t` instead of :math:`k`,
since the truncated level is different from the number of components, which is infinity,
but will be truncated to :math:`T` later.)

.. math::
    \begin{array}{rl}
    q(\mathbf{Z}, \mathbf{V}, \boldsymbol{\mu}, \boldsymbol{\Lambda}) & = q(\mathbf{Z})q(\mathbf{V}, \boldsymbol{\mu}, \boldsymbol{\Lambda}) \\
    q(\mathbf{V}, \boldsymbol{\mu}, \boldsymbol{\Lambda}) & = \prod_{t=1}^{T-1} q(V_t) \prod_{t=1}^T q(\boldsymbol{\mu}_t, \boldsymbol{\Lambda}_t)
    \end{array}

The second term :math:`p(\mathbf{Z}|\mathbf{V})` is slightly different from the first equation
in Blei's paper, page 129, just in terms of representation,
because we still use 1-of-K encoding for :math:`\mathbf{Z}`.
In the truncated variational distribution, we fix a value :math:`T` and :math:`q(V_T = 1) = 1`.
It means that in the :math:`T`-th breaking, we just take the whole of the rest stick. Therefore,

.. math::
    \begin{array}{rl}
    q(\pi_t = 0) & = 1, \quad t > T \\
    q(z_{nt} = 1) & = 0, \quad t > T
    \end{array}

The most part of derivations are exactly the same as in :class:`BayesianGaussianMixture`,
except that all terms involving :math:`\boldsymbol{\pi}`.

In E-step, the updating functions :math:`\mathbf{Z}` are all the same,
except the term

.. math::
    \ln \widetilde{\pi}_k = \mathbb{E}[\ln \pi_k]
    = {}& \mathbb{E}[\ln V_k] + \sum_{i=1}^{k-1} \mathbb{E}[\ln (1- V_i)] \\
    = {}& \psi(\gamma_{k1}) - \psi(\gamma_{k1} + \gamma_{k2}) + \sum_{i=1}^{k-1}
        \left[\psi(\gamma_{i2}) - \psi(\gamma_{i1} + \gamma_{i2})
        \right]

In M-step, we have the updating functions

.. math::
    q^*(V_k) & = \textrm{Beta}(V_k | \gamma_{k1}, \gamma_{k2}) \\
    \gamma_{k1} & = 1 + \sum_{n=1}^N r_{nk} \\
    \gamma_{k2} & = \gamma_0 + \sum_{n=1}^N\sum_{i=k+1}^T r_{nk}


Lower Bound(ELBO)
.................

The lower bound could be estimated by replacing some terms in the lower bound of GMM solved by VB.
In the lower bound of :class:`BayesianGaussianMixture` we just use new :math:`\ln \widetilde{\pi}_k` of DPGMM.

.. math::
    \mathbb{E}[\ln p(\mathbf{Z}|\mathbf{V})] = \sum_{n=1}^N \sum_{k=1}^K r_{nk} \ln \widetilde{\pi}_k

And other new terms we should compute instead of :math:`\mathbb{E} [\ln p(\boldsymbol{\pi})]`
and :math:`\mathbb{E} [\ln q(\boldsymbol{\pi})]` are

.. math::
    \begin{array}{rl}
        \mathbb{E}[\ln p(\mathbf{V})] & = \sum_{t=1}^T \mathbb{E}[\ln p(V_i)]  \\
        & = T \ln F(1, \gamma_0) + (\gamma_0-1)\sum_{t=1}^T[\psi(\gamma_{t2}) - \psi(\gamma_{t1}+ \gamma_{t2})]
    \end{array}

.. math::
  \mathbb{E}[\ln q(\mathbf{V})]
    = \sum_{t=1}^T \left[
                          \ln F(\gamma_{t1}, \gamma_{t2})
                          + (\gamma_{t1}-1) [\psi(\gamma_{t1}) - \psi(\gamma_{t1}+ \gamma_{t2})]
                          + (\gamma_{t2}-1) [\psi(\gamma_{t2}) - \psi(\gamma_{t1}+ \gamma_{t2})]
                      \right]