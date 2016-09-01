.. _bayesain_gaussian_mixture_derivation:

.. currentmodule:: sklearn.mixture

Bayesian Gaussian Mixture Models
================================

The API is identical to that of the :class:`GaussianMixture` class.

The inference algorithm is the one from the following book:

    * `Pattern recognition and machine learning
      <http://www.springer.com/kr/book/9780387310732>`_
      Bishop, Christopher M. Springer, Vol. 4 No. 4, 2006.

While the book presents the parts of the variationnal inference algorithm, it
does not go into detail in the mixture modeling part, which can be just as
complex, or even more. For this reason we present here a full derivation of the
inference algorithm and all the update and lower-bound equations. If you're not
interested in learning how to derive similar algorithms yourself and you're not
interested in changing/debugging the implementation in the scikit this document
is not for you.

This implementation is expected to scale at least as well as EM for
the Gaussian mixture.

Update rules for VB inference
==============================

Here, we present the full mathematical derivation of the update rules
for the :class:`BayesianGaussianMixture` class.

Notation
--------
During this document, we will use the following notation :

- :math:`\mathcal{N}` is the `Normal distribution. <https://en.wikipedia.org/wiki/Normal_distribution>`_
- :math:`\mathcal{W}` is the `Wishart distribution. <https://en.wikipedia.org/wiki/Wishart_distribution>`_
- :math:`Dir` is the `Dirichlet distribution. <https://en.wikipedia.org/wiki/Dirichlet_distribution>`_
- :math:`z_{nk}` is the indicator variable which represents that :math:`\mathbf{x}_n` belongs to the :math:`k`-th components

For convenience, we introduce the following statictics:

.. math::
   :nowrap:

   \begin{eqnarray*}
      N_k & = & \sum_{n=1}^N {r_{nk}} \\
      \bar{\mathbf{x}}_k & = & \frac{1}{N_k} \sum_{n=1}^N {r_{nk} \mathbf{x}_n} \\
      \mathbf{S}_k & = & \frac{1}{N_k} \sum_{n=1}^N {r_{nk} (\mathbf{x}_n - \bar{\mathbf{x}}_k) (\mathbf{x}_n - \bar{\mathbf{x}}_k)^\top}
   \end{eqnarray*}

where :math:`r_{nk}` is the normalized posterior probability that :math:`k` was responsible
for generating :math:`\mathbf{x}_n` (also called responsibility) defined by the
Equation :eq:`vbgmm_responsabilities`.

The full model
--------------
The model of the Bayesian Gaussian Mixture for full covariance matrices is
defined by the joint probability :

.. math::
   :nowrap:
   :label: vbgmm_joint_probability

   \begin{eqnarray*}
      p(\mathbf{X}, \mathbf{Z}, \boldsymbol{\pi}, \boldsymbol{\mu}, \boldsymbol{\Lambda}) & = &
          p(\mathbf{X} | \mathbf{Z}, \boldsymbol{\mu}, \boldsymbol{\Lambda})
          p(\mathbf{Z} | \boldsymbol{\pi})
          p(\boldsymbol{\pi})
          p(\boldsymbol{\mu}|\boldsymbol{\Lambda})
          p(\boldsymbol{\Lambda}) \\
   \end{eqnarray*}

where

.. math::
   :nowrap:

   \begin{eqnarray*}
      p(\mathbf{X} | \mathbf{Z}, \boldsymbol{\mu}, \boldsymbol{\Lambda})
        & = & \prod_{n=1}^N \prod_{k=1}^K \mathcal{N} (\mathbf{x}_n | \boldsymbol{\mu}_k, \boldsymbol{\Lambda}_k^{-1})^{z_{nk}} \\
      p(\mathbf{Z} | \boldsymbol{\pi})
        & = & \prod_{n=1}^N \prod_{k=1}^K \pi_k^{z_{nk}} \\
      p(\boldsymbol{\pi})
        & = & Dir(\boldsymbol{\pi}| \boldsymbol{\alpha}_0) = C(\boldsymbol{\alpha}_0) \prod_{k=1}^K \pi_k^{\alpha_0-1} \\
      p(\boldsymbol{\mu}, \boldsymbol{\Lambda})
        & = & p(\boldsymbol{\mu}|\boldsymbol{\Lambda})p(\boldsymbol{\Lambda})  \\
        & = & \prod_{k=1}^K \mathcal{N}(\boldsymbol{\mu}_k|\mathbf{m}_0, (\beta_0\boldsymbol{\Lambda}_k)^{-1})
                          \mathcal{W}(\boldsymbol{\Lambda}_k| \mathbf{W}_0, \nu_0)
   \end{eqnarray*}

We will consider the variational distribution which factorizes between the
latent variables and the parameters

.. math::
   :nowrap:

   \begin{eqnarray*}
      q(\mathbf{Z}, \boldsymbol{\pi}, \boldsymbol{\mu}, \boldsymbol{\Lambda}) & = &
      q(\mathbf{Z})q(\boldsymbol{\pi}, \boldsymbol{\mu}, \boldsymbol{\Lambda}) \\
      q(\boldsymbol{\pi}, \boldsymbol{\mu}, \boldsymbol{\Lambda}) & = &
      q(\boldsymbol{\pi}) \prod_{k=1}^K q(\boldsymbol{\mu}_k, \boldsymbol{\Lambda}_k)
   \end{eqnarray*}

M-Step
......

Estimate the weights
~~~~~~~~~~~~~~~~~~~~

The goal is to find out the parameters of :math:`\mathcal{N}(\boldsymbol{\mu_k}
| \mathbf{m}_k, (\beta_k \boldsymbol{\Lambda}_k)^{-1})
\mathcal{W}(\boldsymbol{\Lambda}_k|\mathbf{W}_k, \nu_k)`.
Let's take the expectation of the terms about
:math:`\boldsymbol{\pi}, \boldsymbol{\mu}, \boldsymbol{\Lambda}` in the joint
probability :eq:`vbgmm_joint_probability` with the variational distribution of
:math:`\mathbf{Z}`.

.. math::
   :nowrap:
   :label: vbgmm_variationnal_distribution

   \begin{eqnarray*}
      \ln q^\star(\boldsymbol{\pi}, \boldsymbol{\mu}, \boldsymbol{\Lambda})
        & = & \ln p(\boldsymbol{\pi}) + \mathbb{E}_{\mathbf{Z}}[\ln p(\mathbf{Z} | \boldsymbol{\pi})] + \\
        &   & \sum_{k=1}^K \ln p(\boldsymbol{\mu}_k, \boldsymbol{\Lambda}_k) +
            \sum_{k=1}^K \sum_{n=1}^N \mathbb{E}[z_{nk}] \ln \mathcal{N}(\mathbf{x}_n|\boldsymbol{\mu}_k, \boldsymbol{\Lambda}_k) +
            \textrm{const}
   \end{eqnarray*}

Group the terms into only involving :math:`\boldsymbol{\pi}` together with only
involving :math:`\boldsymbol{\mu}` and :math:`\boldsymbol{\Lambda}`, then
variational distribution could be factorized into

.. math::  q(\boldsymbol{\pi}, \boldsymbol{\mu}, \boldsymbol{\Lambda}) =
    q(\boldsymbol{\pi})q(\boldsymbol{\mu}, \boldsymbol{\Lambda})

Taking the term only about :math:`\boldsymbol{\pi}` in the Equation
:eq:`vbgmm_variationnal_distribution`, we obtain

.. math::
   :nowrap:

   \begin{eqnarray*}
      \ln q^\star(\boldsymbol{\pi}) & = & \mathbb{E}_{\mathbf{Z}} [\ln p(\boldsymbol{\pi})]  +
                                      \mathbb{E}_{\mathbf{Z}} [\ln p(\mathbf{Z} | \boldsymbol{\pi})] + \textrm{const} \\
                                & = & (\alpha_0-1)\sum_{k=1}^K \ln \pi_k +
                                      \sum_{k=1}^K\sum_{n=1}^N r_{nk} \ln \pi_k + \textrm{const}
   \end{eqnarray*}

Finally, the weights can be expressed by

.. math::
   :nowrap:
   :label: vbgmm_weights

   \begin{eqnarray*}
      q^\star(\boldsymbol{\pi}) & = & Dir(\boldsymbol{\pi}|\boldsymbol{\alpha}) \\
      \alpha_k & = & \alpha_0 + N_k
   \end{eqnarray*}


Estimate the means
~~~~~~~~~~~~~~~~~~

By taking the term only about :math:`\boldsymbol{\mu}_k` and
:math:`\boldsymbol{\Lambda}_k` in Equation
:eq:`vbgmm_variationnal_distribution`.

.. math::
   :nowrap:

   \begin{eqnarray*}
      \ln q^\star(\boldsymbol{\mu}_k, \boldsymbol{\Lambda}_k)
        & = & \ln \mathcal{N} \left(\boldsymbol{\mu}_k | \mathbf{m}_0,  (\beta_0\boldsymbol{\Lambda}_k)^{-1} \right)
              + \ln \mathcal{W}(\boldsymbol{\Lambda}_k|\mathbf{W}_0, \nu_0) \\
        &   & + \sum_{n=1}^N \mathbb{E}[z_{nk}]\ln \mathcal{N} \left(\mathbf{x}_n|\boldsymbol{\mu}_k, \boldsymbol{\Lambda}_k^{-1} \right)    + \textrm{const} \\
        & = & \frac{1}{2} \ln | \boldsymbol{\Lambda}_k|
              -\frac{\beta_0}{2} (\boldsymbol{\mu}_k - \mathbf{m}_0)^\top \boldsymbol{\Lambda}_k(\boldsymbol{\mu}_k - \mathbf{m}_0) + \\
        &   & \frac{\nu_0 - D - 1}{2} \ln |\boldsymbol{\Lambda}_k| -\frac{1}{2} Tr(\mathbf{W}_0^{-1}\boldsymbol{\Lambda}_k) +   \\
        &   & \frac{1}{2} \left( \sum_{n=1}^N r_{nk} \right) \ln | \boldsymbol{\Lambda}_k |
              - \frac{1}{2}\sum_{n=1}^N r_{nk}
                (\mathbf{x}_n - \boldsymbol{\mu}_k)^\top \boldsymbol{\Lambda}_k (\mathbf{x}_n - \boldsymbol{\mu}_k)
              + \textrm{const}
   \end{eqnarray*}

Using :math:`\ln q^\star(\boldsymbol{\mu}_k, \boldsymbol{\Lambda}_k) =
\ln q^\star(\boldsymbol{\mu}_k | \boldsymbol{\Lambda}_k) +
\ln q^\star(\boldsymbol{\Lambda}_k)`, we first group the terms about
:math:`\boldsymbol{\mu}_k`,

.. math::
   :nowrap:

   \begin{eqnarray*}
      \ln q^\star(\boldsymbol{\mu}_k | \boldsymbol{\Lambda}_k)
       & = & -\frac{1}{2} \boldsymbol{\mu}_k^\top
             \left[\beta_0 + \sum_{n=1}^N r_{nk}  \right] \boldsymbol{\Lambda}_k \boldsymbol{\mu}_k +
             \boldsymbol{\mu}_k^\top \boldsymbol{\Lambda}_k \left[\beta_0 \mathbf{m}_0 + \sum_{n=1}^N r_{nk} \mathbf{x}_n \right] +
             \textrm{const} \\
       & = & -\frac{1}{2} \boldsymbol{\mu}_k^\top [ \beta_0 + N_k ] \boldsymbol{\Lambda}_k \boldsymbol{\mu}_k +
         \boldsymbol{\mu}_k^\top \boldsymbol{\Lambda}_k \left[\beta_0 \mathbf{m}_0 + N_k \bar{\mathbf{x}}_k \right] + \textrm{const}
   \end{eqnarray*}

:math:`q^\star(\boldsymbol{\mu}_k | \boldsymbol{\Lambda}_k)` has a Gaussian
distribution form,

.. math::
   :nowrap:
   :label: vbgmm_means

   \begin{eqnarray*}
      q^\star(\boldsymbol{\mu}_k| \boldsymbol{\Lambda}_k) & = & \mathcal{N} (\boldsymbol{\mu}_k | \mathbf{m}_k, (\beta_k \boldsymbol{\Lambda}_k)^{-1}) \\
      \beta_k & = & \beta_0 + N_k \\
      \mathbf{m}_k & = & \frac{1}{\beta_k}(\beta_0 \mathbf{m}_0 + N_k \bar{\mathbf{x}}_k)
   \end{eqnarray*}


Estimate covariances
~~~~~~~~~~~~~~~~~~~~

Let's estimate :math:`\ln q^\star(\boldsymbol{\Lambda}_k)`

.. math::
   :nowrap:

   \begin{eqnarray*}
      \ln q^\star(\boldsymbol{\Lambda}_k)
         & = & \ln q^\star(\boldsymbol{\mu}_k, \boldsymbol{\Lambda}_k)
               - \ln q^\star(\boldsymbol{\mu}_k | \boldsymbol{\Lambda}_k) \\
         & = & \frac{1}{2} \ln | \boldsymbol{\Lambda}_k|
             - \frac{\beta_0}{2} (\boldsymbol{\mu}_k - \mathbf{m}_0)^\top \boldsymbol{\Lambda}_k(\boldsymbol{\mu}_k - \mathbf{m}_0) +  \\
         &   & \frac{\nu_0 - D - 1}{2} \ln |\boldsymbol{\Lambda}_k|
             - \frac{1}{2} Tr(\mathbf{W}_0^{-1}\boldsymbol{\Lambda}_k) +   \\
         &   & \frac{1}{2} \left( \sum_{n=1}^N r_{nk} \right) \ln | \boldsymbol{\Lambda}_k |
             - \frac{1}{2}\sum_{n=1}^N r_{nk}
                                     (\mathbf{x}_n - \boldsymbol{\mu}_k)^\top \boldsymbol{\Lambda}_k (\mathbf{x}_n - \boldsymbol{\mu}_k) + \\
         &   & -\frac{1}{2}\ln|\boldsymbol{\Lambda}_k| +
             \frac{\beta_k}{2}(\boldsymbol{\mu}_k - \mathbf{m}_k)^\top \boldsymbol{\Lambda}_k (\boldsymbol{\mu}_k - \mathbf{m}_k)
             + \textrm{const} \\
         & = & \frac{\nu_0 + N_k - D - 1}{2}\ln |\boldsymbol{\Lambda}_k| \\
         &   & - \frac{\beta_0}{2} (\boldsymbol{\mu}_k - \mathbf{m}_0)^\top \boldsymbol{\Lambda}_k(\boldsymbol{\mu}_k - \mathbf{m}_0)
            - \frac{1}{2} Tr(\mathbf{W}_0^{-1}\boldsymbol{\Lambda}_k) + \\
         &   & \frac{\beta_k}{2}(\boldsymbol{\mu}_k - \mathbf{m}_k)^\top \boldsymbol{\Lambda}_k(\boldsymbol{\mu}_k - \mathbf{m}_k)
            - \frac{1}{2}\sum_{n=1}^N r_{nk}
                (\mathbf{x}_n - \boldsymbol{\mu}_k)^\top \boldsymbol{\Lambda}_k (\mathbf{x}_n - \boldsymbol{\mu}_k)
            + \textrm{const}
   \end{eqnarray*}

By comparing :math:`\ln | \boldsymbol{\Lambda}_k|$, $Tr(\boldsymbol{\Lambda}_k)`
with :math:`\frac{\nu_k-D-1}{2} \ln |\boldsymbol{\Lambda}_k| - \frac{1}{2}
Tr(\boldsymbol{\Lambda}_k \mathbf{W}_k^{-1})` and applying
:math:`Tr(\mathbf{A}\mathbf{B}) = Tr(\mathbf{B}\mathbf{A})`, we obtain

.. math::
   :nowrap:

   \begin{eqnarray*}
     \nu_k & = & \nu_0 + N_k \\
     \mathbf{W}_k^{-1} & = & \mathbf{W}_0^{-1} +
       \sum_{n=1}^N r_{nk} (\mathbf{x}_n - \boldsymbol{\mu}_k)(\mathbf{x}_n - \boldsymbol{\mu}_k)^\top + \\
      &   & \beta_0 (\boldsymbol{\mu}_k - \mathbf{m}_0)(\boldsymbol{\mu}_k - \mathbf{m}_0)^\top
          - \beta_k (\boldsymbol{\mu}_k - \mathbf{m}_k)(\boldsymbol{\mu}_k - \mathbf{m}_k)^\top
   \end{eqnarray*}

Using the fact that :math:`\sum_{n=1}^N \mathbb{E}[z_{nk}]\mathbf{x}_n
\mathbf{x}_n^\top = N_k \mathbf{S}_k + N_k \bar{\mathbf{x}}_k
\bar{\mathbf{x}}_k^\top`, we have

.. math::
   :nowrap:

   \begin{eqnarray*}
      \mathbf{W}_k^{-1} & = &  \mathbf{W}_0^{-1} + N_k \mathbf{S}_k +
        N_k \bar{\mathbf{x}}_k \bar{\mathbf{x}}_k^\top -
        2 \sum_{n=1}^N r_{nk} \mathbf{x}_n \boldsymbol{\mu} _k^\top
        + N_k \boldsymbol{\mu}_k \boldsymbol{\mu}_k^\top + \\
                        &   & \beta_0 (\boldsymbol{\mu}_k -
        \mathbf{m}_0)(\boldsymbol{\mu}_k - \mathbf{m}_0)^\top -
        \beta_k (\boldsymbol{\mu}_k - \mathbf{m}_k)(\boldsymbol{\mu}_k -
        \mathbf{m}_k)^\top \\
                        & = & \mathbf{W}_0^{-1} + N_k \mathbf{S}_k +
        N_k \bar{\mathbf{x}}_k \bar{\mathbf{x}}_k^\top -
        2 N_k \bar{\mathbf{x}}_k \boldsymbol{\mu}_k^\top +
        N_k \boldsymbol{\mu}_k \boldsymbol{\mu}_k^\top +  \\
                        &   & \beta_0 (\boldsymbol{\mu}_k -
        \mathbf{m}_0)(\boldsymbol{\mu}_k - \mathbf{m}_0)^\top -
        \beta_k (\boldsymbol{\mu}_k - \mathbf{m}_k)(\boldsymbol{\mu}_k -
        \mathbf{m}_k)^\top
   \end{eqnarray*}

Sort out all the terms involing :math:`\boldsymbol{\mu}_k`, the terms
:math:`\bar{\mathbf{x}}_k \boldsymbol{\mu}_k^\top`, :math:`\mathbf{m}_0
\boldsymbol{\mu}_k^\top, \boldsymbol{\mu}_k \boldsymbol{\mu}_k^\top` are all
canceled out using Equations :eq:`vbgmm_means`, it leaves

.. math::
   :nowrap:

   \begin{eqnarray*}
          \mathbf{W}_k^{-1}  & = & \mathbf{W}_0^{-1} + N_k \mathbf{S}_k +
            N_k \bar{\mathbf{x}}_k \bar{\mathbf{x}}_k^\top +
            \beta_0 \mathbf{m}_0 \mathbf{m}_0^\top  -
            \beta_k \mathbf{m}_k \mathbf{m}_k^\top - \\
                             &   & \frac{1}{\beta_k} (\beta_0 +
            \mathbf{m}_0 + N_k \bar{\mathbf{x}}_k)(\beta_0 + \mathbf{m}_0 +
            N_k \bar{\mathbf{x}}_k)^\top \\
                             & = & \mathbf{W}_0^{-1} + N_k \mathbf{S}_k +
            \left(N_k - \frac{N_k^2}{\beta_k} \right) \bar{\mathbf{x}}_k
            \bar{\mathbf{x}}_k^\top - 2 \frac{\beta_0N_k}{\beta_k}
            \mathbf{m}_0 \bar{\mathbf{x}}_k^\top + \left( \beta_0 -
            \frac{\beta_0^2}{\beta_k} \right) \mathbf{m}_0 \mathbf{m}_0^\top
   \end{eqnarray*}

To summarize, we have

.. math::
   :nowrap:
   :label: vbgmm_precisions

   \begin{eqnarray*}
      q^\star(\boldsymbol{\Lambda}_k) & = & \mathcal{W}(\boldsymbol{\Lambda}_k |
        \mathbf{W}_k,  \nu_k) \\
      \nu_k & = & \nu_0 + N_k \\
      \mathbf{W}_k^{-1} & = & \mathbf{W}_0^{-1} + N_k \mathbf{S}_k +
        \frac{N_k \beta_0}{\beta_k}(\bar{\mathbf{x}}_k -
        \mathbf{m}_0)(\bar{\mathbf{x}}_k - \mathbf{m}_0)^\top
   \end{eqnarray*}

Consequently, the M-step is defined by the Equations :eq:`vbgmm_weights`,
:eq:`vbgmm_means` and :eq:`vbgmm_precisions`.


E-Step
......

Estimate Z
~~~~~~~~~~

Using the general result of variational inference and taking the expectation of
the terms about :math:`\mathbf{Z}` in the joint probability with the variational
distribution of :math:`\boldsymbol{\pi}, \boldsymbol{\mu}, \boldsymbol{\Lambda}`,
we have

.. math::
   :nowrap:

   \begin{eqnarray*}
      \ln q^\star(\mathbf{Z}) & = & \mathbb{E}_{\boldsymbol{\pi}}[
        \ln p(\mathbf{Z}|\boldsymbol{\pi})] +
        \mathbb{E}_{\boldsymbol{\mu}, \boldsymbol{\Lambda}}[
        \ln p(\mathbf{X}|\mathbf{Z}, \boldsymbol{\mu}, \boldsymbol{\Lambda})] +
        \textrm{const} \\
                              & = & \sum_{n=1}^N \sum_{k=1}^K z_{nk}
        \ln \rho_{nk} + \textrm{const}
   \end{eqnarray*}

.. math::
   :nowrap:
   :label: vbgmm_log_likelihood

   \begin{eqnarray*}
      \ln \rho_{nk} = \mathbb{E}[\ln \pi_k] +
        \frac{1}{2} \mathbb{E}[\ln |\boldsymbol{\Lambda}_k|] -
          \frac{D}{2} \ln(2\pi) -  \frac{1}{2} \mathbb{E}_{\boldsymbol{\mu}_k,
          \boldsymbol{\Lambda}_k} \left[
            (\mathbf{x}_n - \boldsymbol{\mu}_k)^\top \boldsymbol{\Lambda}_k
            (\mathbf{x}_n - \boldsymbol{\mu}_k)
          \right]
   \end{eqnarray*}

To complete Equation :eq:`vbgmm_log_likelihood`, we need to compute three
expectations.

By using the property of covariance matrix :math:`(\beta_k
\boldsymbol{\Lambda}_k)^{-1} =
\mathbb{E}[\boldsymbol{\mu} \boldsymbol{\mu}^\top] -
\mathbb{E}[\boldsymbol{\mu}] \mathbb{E}[\boldsymbol{\mu}]^\top` we obtain

.. math::
   :nowrap:

   \begin{eqnarray*}
     \mathbb{E}_{\boldsymbol{\mu}_k, \boldsymbol{\Lambda}_k}
         \left[
           (\mathbf{x}_n - \boldsymbol{\mu}_k)^\top \boldsymbol{\Lambda}_k
           (\mathbf{x}_n - \boldsymbol{\mu}_k)
         \right]
       & = & \mathbb{E}_{\boldsymbol{\Lambda}_k}
         \left[ \mathbb{E}_{\boldsymbol{\mu}_k}
             \left[
                 Tr \left[
                       \boldsymbol{\Lambda}_k
                       (\mathbf{x}_n\mathbf{x}_n^\top
                       - \mathbf{x}_n \boldsymbol{\mu}_k^\top
                       - \boldsymbol{\mu}_k \mathbf{x}_n^\top
                       + \boldsymbol{\mu}_k \boldsymbol{\mu}_k^\top)
                    \right]
             \right]
         \right] \\
       & = & \mathbb{E}_{\boldsymbol{\Lambda}_k}
         \left[
            Tr \left[
                  \boldsymbol{\Lambda}_k (\mathbf{x}_n\mathbf{x}_n^\top -
                  \mathbf{x}_n \mathbf{m}_k^\top -
                  \mathbf{m}_k \mathbf{x}_n^\top +
                  (\beta_k \boldsymbol{\Lambda}_k)^{-1} +
                  \mathbf{m}_k \mathbf{m}_k^\top )
              \right]
         \right] \\
       & = & D \beta_k^{-1} + \nu_k (\mathbf{x}_n -
         \mathbf{m}_k)^\top \mathbf{W}_k (\mathbf{x}_n - \mathbf{m}_k)
   \end{eqnarray*}

The three expectations are finally defined by :

.. math::
   :nowrap:
   :label: vbgmm_expectation

   \begin{eqnarray*}
      \ln \widetilde{\pi}_k & = & \mathbb{E}[\ln \pi_k] = \psi(\alpha_k) -
        \psi(\hat{\alpha}) \\
      \ln \widetilde{\Lambda}_k & = &
        \mathbb{E}[\ln |\boldsymbol{\Lambda}_k|] =
          \sum_{i=1}^D \psi \left( \frac{\nu_k + 1 - i}{2} \right) +
          D \ln 2 + \ln |\mathbf{W}_k| \\
      \mathbb{E}_{\boldsymbol{\mu}_k, \boldsymbol{\Lambda}_k}
        \left[
          (\mathbf{x}_n - \boldsymbol{\mu}_k)^\top \boldsymbol{\Lambda}_k
          (\mathbf{x}_n - \boldsymbol{\mu}_k)
        \right]
      & = & D \beta_k^{-1} + \nu_k (\mathbf{x}_n -
        \mathbf{m}_k)^\top \mathbf{W}_k (\mathbf{x}_n - \mathbf{m}_k)
   \end{eqnarray*}

After normalizing :math:`\rho_{nk}`, we have responsibility for each :math:`n`
for each :math:`k` and the posterior multinomial distribution.

.. math::
   :nowrap:
   :label: vbgmm_responsabilities

   \begin{eqnarray*}
      q^\star(\mathbf{Z}) & = & \prod_{n=1}^N\prod_{k=1}^K r_{nk}^{z_{nk}} \\
      r_{nk} & = & \frac{\rho_{nk}}{\sum_{i=1}^K \rho_{ni}}
   \end{eqnarray*}

Lower bound
...........
The lower bound :math:`\mathcal{L}` is defined by

.. math::
   :nowrap:
   :label: vbgmm_responsabilities

   \begin{eqnarray*}
     \mathcal{L}
       & = &  \sum_{\mathbf{Z}} \iiint q(\mathbf{Z}, \boldsymbol{\pi},
         \boldsymbol{\mu}, \boldsymbol{\Lambda})
         \ln \left\{
           \frac{p(\mathbf{X}, \mathbf{Z}, \boldsymbol{\pi}, \boldsymbol{\mu},
             \boldsymbol{\Lambda})} {q(\mathbf{Z}, \boldsymbol{\pi},
             \boldsymbol{\mu}, \boldsymbol{\Lambda})}\right\}
           \mathrm{d} \boldsymbol{\pi} \mathrm{d} \boldsymbol{\mu}
           \mathrm{d} \boldsymbol{\Lambda} \\
       & = &  \mathbb{E}[\ln p(\mathbf{X}|\mathbf{Z},\boldsymbol{\mu},
         \boldsymbol{\Lambda})] +
         \mathbb{E} [\ln p(\mathbf{Z}|\boldsymbol{\pi})] +
         \mathbb{E} [\ln p(\boldsymbol{\pi})] +
         \mathbb{E} [\ln p(\boldsymbol{\mu}, \boldsymbol{\Lambda})] \\
       &   &  -\mathbb{E} [\ln q(\mathbf{Z})] -
         \mathbb{E} [\ln q(\boldsymbol{\pi})] -
         \mathbb{E} [\ln q(\boldsymbol{\mu}, \boldsymbol{\Lambda})]
   \end{eqnarray*}

Estimate lower bound terms
--------------------------
To compute :math:`\mathbb{E}[\ln p(\mathbf{X}|\mathbf{Z},\boldsymbol{\mu},
\boldsymbol{\Lambda})]` we can use the Equation :eq:`vbgmm_expectation` and use
the fact :math:`\sum_{n=1}^N \mathbb{E}[z_{nk}] \mathbf{x}_n \mathbf{x}_n^\top =
N_k \mathbf{S}_k + N_k \bar{\mathbf{x}}_k \bar{\mathbf{x}}_k^\top`

.. math::
   :nowrap:

   \begin{eqnarray*}
      \mathbb{E}[\ln p(\mathbf{X}|\mathbf{Z},\boldsymbol{\mu},
        \boldsymbol{\Lambda})]
      & = & \frac{1}{2} \mathbb{E}
            \left[
              \sum_{n=1}^N \sum_{k=1}^K z_{nk}
                \left\{
                  \ln |\boldsymbol{\Lambda}_k| - D \ln (2\pi)
                  - (\mathbf{x}_n - \boldsymbol{\mu}_k)^\top
                  \boldsymbol{\Lambda}_k (\mathbf{x}_n - \boldsymbol{\mu}_k)
                \right\}
            \right] \\
      & = & \frac{1}{2} \sum_{k=1}^K N_k \left\{ \ln \widetilde{\Lambda}_k -
        D \ln(2\pi) \right\} - \sum_{n=1}^N \sum_{k=1}^K r_{nk} [
          D\beta_k^{-1} + \nu_k (\mathbf{x}_n - \mathbf{m}_k)^\top \mathbf{W}_k
          (\mathbf{x}_n - \mathbf{m}_k)] \\
      & = & \frac{1}{2} \sum_{k=1}^K N_k
        \left\{
          \ln \widetilde{\Lambda}_k - D \ln(2\pi) - D\beta_k^{-1} -
          \nu_k Tr[\mathbf{S}_k \mathbf{W}_k] - \nu_k (\bar{\mathbf{x}}_k -
          \mathbf{m}_k)^\top \mathbf{W}_k (\bar{\mathbf{x}}_k - \mathbf{m}_k)]
        \right\}
   \end{eqnarray*}

The term :math:`\mathbb{E} [\ln p(\boldsymbol{\mu}, \boldsymbol{\Lambda})]`
is defined from the Equation :eq:`vbgmm_expectation` by

.. math::
   :nowrap:

   \begin{eqnarray*}
      \mathbb{E} [\ln p(\boldsymbol{\mu}, \boldsymbol{\Lambda})]
      & = & \sum_{k=1}^K \mathbb{E}
        \left[
          \ln \mathcal{N}(\boldsymbol{\mu}_k |
            \mathbf{m}_0, (\beta_0 \boldsymbol{\Lambda}_k)^{-1}) +
          \ln \mathcal{W}(\boldsymbol{\Lambda}_k | \mathbf{W}_0, \nu_0 )
        \right] \\
      & = & \frac{1}{2} \sum_{k=1}^K
          \left\{
              D \ln \frac{\beta_0}{2\pi}
              + \ln \widetilde{\Lambda}_k
              - \frac{D \beta_0}{\beta_k}
              - \beta_0 \nu_k (\mathbf{m}_k - \mathbf{m}_0)^\top
                \mathbf{W}_k (\mathbf{m}_k - \mathbf{m}_0)
          \right\} \\
      &   & + K \ln B(\mathbf{W}_0, \nu_0) +
        \frac{\nu_0-D-1}{2} \sum_{k=1}^K \ln \widetilde{\Lambda}_k -
        \frac{1}{2} \sum_{k=1}^K \nu_k Tr(\mathbf{W}_0^{-1}\mathbf{W}_k)
   \end{eqnarray*}

We can estimate
:math:`\mathbb{E}[\ln q(\boldsymbol{\mu}, \boldsymbol{\Lambda})]` by

.. math::
   :nowrap:

   \begin{eqnarray*}
     \mathbb{E}[\ln q(\boldsymbol{\mu}, \boldsymbol{\Lambda})]
     & = & \sum_{k=1}^K \mathbb{E}[
       \ln \mathcal{N}(\boldsymbol{\mu}_k | \mathbf{m}_k,
         (\beta_k \boldsymbol{\Lambda}_k)^{-1}) +
       \ln \mathcal{W}(\boldsymbol{\Lambda}_k | \mathbf{W}_k, \nu_k) ] \\
     & = & \sum_{k=1}^K
       \left\{
         \frac{1}{2} \ln \widetilde{\Lambda}_k
         + \frac{D}{2} \ln \frac{\beta_k}{2 \pi}
         - \frac{D}{2}
         -\ln B(\mathbf{W}, \nu)
         - \frac{(\nu - D - 1)}{2} \mathbb{E}[\ln |\boldsymbol{\Lambda}|]
         + \frac{\nu D}{2}]
       \right\}
   \end{eqnarray*}

Lower bound combination and simplifications
...........................................

Some simplifications and combination of terms can be performed in the Equation
:eq:`vbgmm_responsabilities`:

.. math::
   :nowrap:

   \begin{eqnarray*}
     \mathbb{E}[\ln p(\mathbf{Z}|\boldsymbol{\pi})]
        & = & \sum_{n=1}^N \sum_{k=1}^K r_{nk} \ln \widetilde{\pi}_k =
          \sum_{k=1}^K N_k \ln \widetilde{\pi}_k \\
     \mathbb{E}[\ln p(\boldsymbol{\pi})]
        & = & \ln C(\boldsymbol{\alpha}_0) +
          (\alpha_0 - 1) \sum_{k=1}^K \ln \widetilde{\pi}_k \\
     \mathbb{E}[\ln q(\boldsymbol{\pi})]
        & = & \sum_{k=1}^K(\alpha_k - 1) \ln \tilde{\pi}_k +
          \ln C(\boldsymbol{\alpha})
   \end{eqnarray*}

Thus,

.. math::
   :nowrap:

   \begin{eqnarray*}
     \mathbb{E}[\ln p(\mathbf{Z}|\boldsymbol{\pi})] +
       \mathbb{E}[\ln p(\boldsymbol{\pi})] -
       \mathbb{E}[\ln q(\boldsymbol{\pi})]
     & = & \ln C(\boldsymbol{\alpha}_0) - \ln C(\boldsymbol{\alpha}) +
       \sum_{k=1}^K(\alpha_0 + N_k - \alpha_k - 1 + 1) \ln \tilde{\pi}_k \\
     & = & \ln C(\boldsymbol{\alpha}_0) - \ln C(\boldsymbol{\alpha})
   \end{eqnarray*}

In a similar way we can simplify

.. math::
   :nowrap:

   \begin{eqnarray*}
     \mathbb{E}[\ln p(\mathbf{X}|\mathbf{Z},\boldsymbol{\mu},
         \boldsymbol{\Lambda})] +
       \mathbb{E}[\ln p(\boldsymbol{\mu}, \boldsymbol{\Lambda})] -
       \mathbb{E}[\ln q(\boldsymbol{\mu}, \boldsymbol{\Lambda})]
     & = & \frac{K}{2} \ln B(\mathbf{W_0}, \nu_0) -
       \frac{1}{2} \sum_{k=1}^K \ln B(\mathbf{W}, \nu) + \\
     &   & K \ln \beta_0 - \sum_{k=1}^K \ln \beta_k - \frac{ND}{2} \ln 2\pi
     \end{eqnarray*}

Consequently the lower value is :

.. math::
   :nowrap:
   :label: vbgmm_responsabilities

   \begin{eqnarray*}
     \mathcal{L} & = & \ln C(\boldsymbol{\alpha}_0) -
                       \ln C(\boldsymbol{\alpha}) +
                       K \ln \beta_0 - \sum_{k=1}^K \ln \beta_k + \\
                 &   & \frac{K}{2} \ln B(\mathbf{W_0}, \nu_0) -
                       \frac{1}{2} \sum_{k=1}^K \ln B(\mathbf{W}, \nu) - \\
                 &   & \sum_{n=1}^N \sum_{k=1}^K r_{nk} \ln r_{nk} -
                       \frac{ND}{2} \ln 2 \pi
   \end{eqnarray*}


References
----------
[1] `Bishop, Christopher M. (2006). "Pattern recognition and machine learning".
     Vol. 4 No. 4. New York: Springer. <http://www.springer.com/kr/book/9780387310732>`_

[2] `Hagai Attias. (2000). "A Variational Bayesian Framework for Graphical Models".
     In Advances in Neural Information Processing Systems 12
     <http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.36.2841&rep=rep1&type=pdf>`_
