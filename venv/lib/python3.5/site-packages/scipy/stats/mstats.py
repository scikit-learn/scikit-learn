"""
===================================================================
Statistical functions for masked arrays (:mod:`scipy.stats.mstats`)
===================================================================

.. currentmodule:: scipy.stats.mstats

This module contains a large number of statistical functions that can
be used with masked arrays.

Most of these functions are similar to those in scipy.stats but might
have small differences in the API or in the algorithm used. Since this
is a relatively new package, some API changes are still possible.

.. autosummary::
   :toctree: generated/

   argstoarray
   chisquare
   count_tied_groups
   describe
   f_oneway
   find_repeats
   friedmanchisquare
   kendalltau
   kendalltau_seasonal
   kruskalwallis
   ks_twosamp
   kurtosis
   kurtosistest
   linregress
   mannwhitneyu
   plotting_positions
   mode
   moment
   mquantiles
   msign
   normaltest
   obrientransform
   pearsonr
   plotting_positions
   pointbiserialr
   rankdata
   scoreatpercentile
   sem
   skew
   skewtest
   spearmanr
   theilslopes
   tmax
   tmean
   tmin
   trim
   trima
   trimboth
   trimmed_stde
   trimr
   trimtail
   tsem
   ttest_onesamp
   ttest_ind
   ttest_onesamp
   ttest_rel
   tvar
   variation
   winsorize
   zmap
   zscore
   compare_medians_ms
   gmean
   hdmedian
   hdquantiles
   hdquantiles_sd
   hmean
   idealfourths
   kruskal
   ks_2samp
   median_cihs
   meppf
   mjci
   mquantiles_cimj
   rsh
   sen_seasonal_slopes
   trimmed_mean
   trimmed_mean_ci
   trimmed_std
   trimmed_var
   ttest_1samp

"""
from __future__ import division, print_function, absolute_import

from .mstats_basic import *
from .mstats_extras import *
# Functions that support masked array input in stats but need to be kept in the
# mstats namespace for backwards compatibility:
from scipy.stats import gmean, hmean, zmap, zscore, chisquare
