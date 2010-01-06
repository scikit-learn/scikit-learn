#! /usr/bin/env python
# Last Change: Mon Feb 02 11:00 PM 2009 J

# Various utilities for examples 

import numpy as N

"""Different tools for pre processing, like whitening or scaling data."""

_DEF_SCALE_MODE = 'sym'

#---------------------------------------------------------------------------
# Definition of scaling-related function (handle each feature independantly)
#---------------------------------------------------------------------------
def _scale_factor(data, mode = _DEF_SCALE_MODE):
    """Compute the scaling factors for data to be normalized.
    
    Note
    ----
    
    Does not handle data with nan."""
    n = N.min(data, 0)
    m = N.max(data, 0)
    if mode == 'sym':
        t = n + 0.5 * (m - n)
        s = 0.5 * (m - n)
    elif mode == 'right':
        t = n
        s = m - n
    else:
        raise ValueError("Mode %s not recognized" % mode)
    return t, s
    
def _nan_scale_factor(data, mode = _DEF_SCALE_MODE):
    """Compute the scaling factors for data to be normalized.
    
    This version handled data with Nan."""
    n = N.nanmin(data, 0)
    m = N.nanmax(data, 0)
    # If any t or s has Nan, this means that one feature has only Nan. This
    # will propagate Nan to new preprocessed data, which does not really
    # make sense: raise an exception in this case.
    if not (N.all(N.isfinite(n)) and N.all(N.isfinite(m))):
        raise ValueError("Nan scale factors: is any feature of data full"\
                         "of Nan ?")
    if mode == 'sym':
        t = n + 0.5 * (m - n)
        s = 0.5 * (m - n)
    elif mode == 'right':
        t = n
        s = m - n
    else:
        raise ValueError("Mode %s not recognized" % mode)
    return t, s
    
def schandlenan(f, note):
    """Decorator to generate scaling function handling / not handling Nan."""
    def decorator(func):
        ""
        def wrapper(data, mode = _DEF_SCALE_MODE):
            ""
            return func(data, mode, f)
        wrapper.__name__ = func.__name__
        wrapper.__dict__ = func.__dict__
        wrapper.__doc__ = \
"""Linearly scale data in place such as each col is in the range
[0..1] or [-1..1].

:Parameters:
    data : ndarray
        the data to scale. One feature per column (eg the
        normalization is done column-wise).
    mode : string
        - 'sym': normalized data are in the range [-1..1]
        - 'right': normalized data are in the range [0..1]

:Returns:
    s : ndarray
        the scaling factor.
    t : ndarray
        the translation factor
        
:SeeAlso:

Scaler, whiten.

Note:
-----

You can retrieve the original values with data = s * scaled + t.  

This is intended to normalize data for pre processing; in
perticular, the range after normalized do not have strong
constraints: some values may be higher than 1 due to precision
problems.\n\n""" + note
        return wrapper
    return decorator

@schandlenan(_scale_factor, "This function does NOT handle Nan.")
def scale(data, mode, scf):
    t, s = scf(data, mode)
    data -= t
    data /= s
    return s, t

@schandlenan(_nan_scale_factor, "This function DOES handle Nan.")
def nanscale(data, mode, scf):
    t, s = scf(data, mode)
    data -= t
    data /= s
    return s, t
#def scale(data, mode = _DEF_SCALE_MODE):
#    """Linearly scale data in place such as each col is in the range [0..1] or
#    [-1..1].
#
#    :Parameters:
#        data : ndarray
#            the data to scale. One feature per column (eg the normalization is
#            done column-wise).
#        mode : string
#            - 'sym': normalized data are in the range [-1..1]
#            - 'right': normalized data are in the range [0..1]
#
#    :Returns:
#        s : ndarray
#            the scaling factor.
#        t : ndarray
#            the translation factor
#            
#    :SeeAlso:
#
#    Scaler, whiten.
#
#    Note:
#    -----
#
#    Handle data with Nan values (are ignored)
#
#    You can retrieve the original values with data = s * scaled + t.  
#
#    This is intended to normalize data for pre processing; in perticular, the
#    range after normalized do not have strong constraints: some values may be
#    higher than 1 due to precision problems."""
#    t, s = _scale_factor(data, mode)
#    data -= t
#    data /= s
#    return s, t

def whiten():
    """Whiten data."""
    raise NotImplementedError("whitening not implemented yet")

class Scaler:
    """Class to implement a scaler, eg an object which can scale data, keep the
    scale factors, and rescale/unscale further data using those factors.

    For example, in supervised training, you may want to rescale some data,
    usually using the training data. Once the training is done with the scaled
    data, you need to scale the testing and unknown data by the same factor,
    and maybe also to "unscale" them afterwards."""
    def __init__(self, data, mode = _DEF_SCALE_MODE):
        """Init the scaler, computing scaling factor from data."""
        t, s = _scale_factor(data, mode)
        self.t = t
        self.s = s

    def scale(self, data):
        """Scale data *in-place*."""
        try:
            data -= self.t
            data /= self.s
        except ValueError, e:
            raise ValueError("data to scale should have the same number of"\
                             "features than data used when initializing the"\
                             "scaler")
        return data

    def unscale(self, data):
        """Unscale data *in-place*."""
        try:
            data *= self.s
            data += self.t
        except ValueError, e:
            raise ValueError("data to unscale should have the same number of"\
                             "features than data used when initializing the"\
                             "scaler")
        return data

class NanScaler(Scaler):
    """Special scaler which ignore Nan data."""
    def __init__(self, data, mode = _DEF_SCALE_MODE):
        """Init the scaler, computing scaling factor from data.

        Note
        ----

        If one feature has only Nan, the scaling parameters will contain Nan
        values, and as such, will propagate Nan to newly data fed to
        preprocess. To avoid this, an exception is raised."""
        t, s = _nan_scale_factor(data, mode)
        self.t = t
        self.s = s
