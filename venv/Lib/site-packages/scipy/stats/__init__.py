"""
.. _statsrefmanual:

==========================================
Statistical functions (:mod:`scipy.stats`)
==========================================

.. currentmodule:: scipy.stats

This module contains a large number of probability distributions as
well as a growing library of statistical functions.

Each univariate distribution is an instance of a subclass of `rv_continuous`
(`rv_discrete` for discrete distributions):

.. autosummary::
   :toctree: generated/

   rv_continuous
   rv_discrete
   rv_histogram

Continuous distributions
========================

.. autosummary::
   :toctree: generated/

   alpha             -- Alpha
   anglit            -- Anglit
   arcsine           -- Arcsine
   argus             -- Argus
   beta              -- Beta
   betaprime         -- Beta Prime
   bradford          -- Bradford
   burr              -- Burr (Type III)
   burr12            -- Burr (Type XII)
   cauchy            -- Cauchy
   chi               -- Chi
   chi2              -- Chi-squared
   cosine            -- Cosine
   crystalball       -- Crystalball
   dgamma            -- Double Gamma
   dweibull          -- Double Weibull
   erlang            -- Erlang
   expon             -- Exponential
   exponnorm         -- Exponentially Modified Normal
   exponweib         -- Exponentiated Weibull
   exponpow          -- Exponential Power
   f                 -- F (Snecdor F)
   fatiguelife       -- Fatigue Life (Birnbaum-Saunders)
   fisk              -- Fisk
   foldcauchy        -- Folded Cauchy
   foldnorm          -- Folded Normal
   genlogistic       -- Generalized Logistic
   gennorm           -- Generalized normal
   genpareto         -- Generalized Pareto
   genexpon          -- Generalized Exponential
   genextreme        -- Generalized Extreme Value
   gausshyper        -- Gauss Hypergeometric
   gamma             -- Gamma
   gengamma          -- Generalized gamma
   genhalflogistic   -- Generalized Half Logistic
   geninvgauss       -- Generalized Inverse Gaussian
   gilbrat           -- Gilbrat
   gompertz          -- Gompertz (Truncated Gumbel)
   gumbel_r          -- Right Sided Gumbel, Log-Weibull, Fisher-Tippett, Extreme Value Type I
   gumbel_l          -- Left Sided Gumbel, etc.
   halfcauchy        -- Half Cauchy
   halflogistic      -- Half Logistic
   halfnorm          -- Half Normal
   halfgennorm       -- Generalized Half Normal
   hypsecant         -- Hyperbolic Secant
   invgamma          -- Inverse Gamma
   invgauss          -- Inverse Gaussian
   invweibull        -- Inverse Weibull
   johnsonsb         -- Johnson SB
   johnsonsu         -- Johnson SU
   kappa4            -- Kappa 4 parameter
   kappa3            -- Kappa 3 parameter
   ksone             -- Distribution of Kolmogorov-Smirnov one-sided test statistic
   kstwo             -- Distribution of Kolmogorov-Smirnov two-sided test statistic
   kstwobign         -- Limiting Distribution of scaled Kolmogorov-Smirnov two-sided test statistic.
   laplace           -- Laplace
   laplace_asymmetric    -- Asymmetric Laplace
   levy              -- Levy
   levy_l
   levy_stable
   logistic          -- Logistic
   loggamma          -- Log-Gamma
   loglaplace        -- Log-Laplace (Log Double Exponential)
   lognorm           -- Log-Normal
   loguniform        -- Log-Uniform
   lomax             -- Lomax (Pareto of the second kind)
   maxwell           -- Maxwell
   mielke            -- Mielke's Beta-Kappa
   moyal             -- Moyal
   nakagami          -- Nakagami
   ncx2              -- Non-central chi-squared
   ncf               -- Non-central F
   nct               -- Non-central Student's T
   norm              -- Normal (Gaussian)
   norminvgauss      -- Normal Inverse Gaussian
   pareto            -- Pareto
   pearson3          -- Pearson type III
   powerlaw          -- Power-function
   powerlognorm      -- Power log normal
   powernorm         -- Power normal
   rdist             -- R-distribution
   rayleigh          -- Rayleigh
   rice              -- Rice
   recipinvgauss     -- Reciprocal Inverse Gaussian
   semicircular      -- Semicircular
   skewnorm          -- Skew normal
   t                 -- Student's T
   trapezoid         -- Trapezoidal
   triang            -- Triangular
   truncexpon        -- Truncated Exponential
   truncnorm         -- Truncated Normal
   tukeylambda       -- Tukey-Lambda
   uniform           -- Uniform
   vonmises          -- Von-Mises (Circular)
   vonmises_line     -- Von-Mises (Line)
   wald              -- Wald
   weibull_min       -- Minimum Weibull (see Frechet)
   weibull_max       -- Maximum Weibull (see Frechet)
   wrapcauchy        -- Wrapped Cauchy

Multivariate distributions
==========================

.. autosummary::
   :toctree: generated/

   multivariate_normal    -- Multivariate normal distribution
   matrix_normal          -- Matrix normal distribution
   dirichlet              -- Dirichlet
   wishart                -- Wishart
   invwishart             -- Inverse Wishart
   multinomial            -- Multinomial distribution
   special_ortho_group    -- SO(N) group
   ortho_group            -- O(N) group
   unitary_group          -- U(N) group
   random_correlation     -- random correlation matrices
   multivariate_t         -- Multivariate t-distribution
   multivariate_hypergeom -- Multivariate hypergeometric distribution

Discrete distributions
======================

.. autosummary::
   :toctree: generated/

   bernoulli         -- Bernoulli
   betabinom         -- Beta-Binomial
   binom             -- Binomial
   boltzmann         -- Boltzmann (Truncated Discrete Exponential)
   dlaplace          -- Discrete Laplacian
   geom              -- Geometric
   hypergeom         -- Hypergeometric
   logser            -- Logarithmic (Log-Series, Series)
   nbinom            -- Negative Binomial
   nhypergeom        -- Negative Hypergeometric
   planck            -- Planck (Discrete Exponential)
   poisson           -- Poisson
   randint           -- Discrete Uniform
   skellam           -- Skellam
   zipf              -- Zipf
   yulesimon         -- Yule-Simon

An overview of statistical functions is given below.
Several of these functions have a similar version in
`scipy.stats.mstats` which work for masked arrays.

Summary statistics
==================

.. autosummary::
   :toctree: generated/

   describe          -- Descriptive statistics
   gmean             -- Geometric mean
   hmean             -- Harmonic mean
   kurtosis          -- Fisher or Pearson kurtosis
   mode              -- Modal value
   moment            -- Central moment
   skew              -- Skewness
   kstat             --
   kstatvar          --
   tmean             -- Truncated arithmetic mean
   tvar              -- Truncated variance
   tmin              --
   tmax              --
   tstd              --
   tsem              --
   variation         -- Coefficient of variation
   find_repeats
   trim_mean
   gstd              -- Geometric Standard Deviation
   iqr
   sem
   bayes_mvs
   mvsdist
   entropy
   median_absolute_deviation
   median_abs_deviation

Frequency statistics
====================

.. autosummary::
   :toctree: generated/

   cumfreq
   itemfreq
   percentileofscore
   scoreatpercentile
   relfreq

.. autosummary::
   :toctree: generated/

   binned_statistic     -- Compute a binned statistic for a set of data.
   binned_statistic_2d  -- Compute a 2-D binned statistic for a set of data.
   binned_statistic_dd  -- Compute a d-D binned statistic for a set of data.

Correlation functions
=====================

.. autosummary::
   :toctree: generated/

   f_oneway
   pearsonr
   spearmanr
   pointbiserialr
   kendalltau
   weightedtau
   linregress
   siegelslopes
   theilslopes
   multiscale_graphcorr

Statistical tests
=================

.. autosummary::
   :toctree: generated/

   ttest_1samp
   ttest_ind
   ttest_ind_from_stats
   ttest_rel
   chisquare
   cramervonmises
   power_divergence
   kstest
   ks_1samp
   ks_2samp
   epps_singleton_2samp
   mannwhitneyu
   tiecorrect
   rankdata
   ranksums
   wilcoxon
   kruskal
   friedmanchisquare
   brunnermunzel
   combine_pvalues
   jarque_bera

.. autosummary::
   :toctree: generated/

   ansari
   bartlett
   levene
   shapiro
   anderson
   anderson_ksamp
   binom_test
   fligner
   median_test
   mood
   skewtest
   kurtosistest
   normaltest

Transformations
===============

.. autosummary::
   :toctree: generated/

   boxcox
   boxcox_normmax
   boxcox_llf
   yeojohnson
   yeojohnson_normmax
   yeojohnson_llf
   obrientransform
   sigmaclip
   trimboth
   trim1
   zmap
   zscore

Statistical distances
=====================

.. autosummary::
   :toctree: generated/

   wasserstein_distance
   energy_distance

Random variate generation
=========================

.. autosummary::
   :toctree: generated/

   rvs_ratio_uniforms

Circular statistical functions
==============================

.. autosummary::
   :toctree: generated/

   circmean
   circvar
   circstd

Contingency table functions
===========================

.. autosummary::
   :toctree: generated/

   chi2_contingency
   contingency.expected_freq
   contingency.margins
   fisher_exact

Plot-tests
==========

.. autosummary::
   :toctree: generated/

   ppcc_max
   ppcc_plot
   probplot
   boxcox_normplot
   yeojohnson_normplot


Masked statistics functions
===========================

.. toctree::

   stats.mstats


Univariate and multivariate kernel density estimation
=====================================================

.. autosummary::
   :toctree: generated/

   gaussian_kde

Warnings used in :mod:`scipy.stats`
===================================

.. autosummary::
   :toctree: generated/

   F_onewayConstantInputWarning
   F_onewayBadInputSizesWarning
   PearsonRConstantInputWarning
   PearsonRNearConstantInputWarning
   SpearmanRConstantInputWarning

For many more stat related functions install the software R and the
interface package rpy.

"""
from .stats import *
from .distributions import *
from .morestats import *
from ._binned_statistic import *
from .kde import gaussian_kde
from . import mstats
from .contingency import chi2_contingency
from ._multivariate import *

__all__ = [s for s in dir() if not s.startswith("_")]  # Remove dunders.

from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
