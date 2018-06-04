# -*- coding: utf-8 -*-
"""Image restoration module.

"""

from .deconvolution import wiener, unsupervised_wiener, richardson_lucy
from .unwrap import unwrap_phase
from ._denoise import (denoise_tv_chambolle, denoise_tv_bregman,
                       denoise_bilateral, denoise_wavelet, estimate_sigma)
from ._cycle_spin import cycle_spin
from .non_local_means import denoise_nl_means
from .inpaint import inpaint_biharmonic


__all__ = ['wiener',
           'unsupervised_wiener',
           'richardson_lucy',
           'unwrap_phase',
           'denoise_tv_bregman',
           'denoise_tv_chambolle',
           'denoise_bilateral',
           'denoise_wavelet',
           'denoise_nl_means',
           'estimate_sigma',
           'inpaint_biharmonic',
           'cycle_spin']
