"""
Test used to verify PyWavelets Discrete Wavelet Transform computation
accuracy against MathWorks Wavelet Toolbox.
"""

from __future__ import division, print_function, absolute_import

import os
import numpy as np
from numpy.testing import assert_, dec, run_module_suite

import pywt

if 'PYWT_XSLOW' in os.environ:
    # Run a more comprehensive set of problem sizes.  This could take more than
    # an hour to complete.
    size_set = 'full'
    use_precomputed = False
else:
    size_set = 'reduced'
    use_precomputed = True

if use_precomputed:
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    matlab_data_file = os.path.join(data_dir, 'dwt_matlabR2012a_result.npz')
    matlab_result_dict = np.load(matlab_data_file)
else:
    try:
        from pymatbridge import Matlab
        mlab = Matlab()
        _matlab_missing = False
    except ImportError:
        print("To run Matlab compatibility tests you need to have MathWorks "
              "MATLAB, MathWorks Wavelet Toolbox and the pymatbridge Python "
              "package installed.")
        _matlab_missing = True

# list of mode names in pywt and matlab
modes = [('zero', 'zpd'),
         ('constant', 'sp0'),
         ('symmetric', 'sym'),
         ('reflect', 'symw'),
         ('periodic', 'ppd'),
         ('smooth', 'sp1'),
         ('periodization', 'per')]

families = ('db', 'sym', 'coif', 'bior', 'rbio')
wavelets = sum([pywt.wavelist(name) for name in families], [])


def _get_data_sizes(w):
    """ Return the sizes to test for wavelet w. """
    if size_set == 'full':
        data_sizes = list(range(w.dec_len, 40)) + \
            [100, 200, 500, 1000, 50000]
    else:
        data_sizes = (w.dec_len, w.dec_len + 1)
    return data_sizes


@dec.skipif(use_precomputed or _matlab_missing)
@dec.slow
def test_accuracy_pymatbridge():
    rstate = np.random.RandomState(1234)
    # max RMSE (was 1.0e-10, is reduced to 5.0e-5 due to different coefficents)
    epsilon = 5.0e-5
    epsilon_pywt_coeffs = 1.0e-10
    mlab.start()
    try:
        for wavelet in wavelets:
            w = pywt.Wavelet(wavelet)
            mlab.set_variable('wavelet', wavelet)
            for N in _get_data_sizes(w):
                data = rstate.randn(N)
                mlab.set_variable('data', data)
                for pmode, mmode in modes:
                    ma, md = _compute_matlab_result(data, wavelet, mmode)
                    yield _check_accuracy, data, w, pmode, ma, md, wavelet, epsilon
                    ma, md = _load_matlab_result_pywt_coeffs(data, wavelet, mmode)
                    yield _check_accuracy, data, w, pmode, ma, md, wavelet, epsilon_pywt_coeffs

    finally:
        mlab.stop()


@dec.skipif(not use_precomputed)
@dec.slow
def test_accuracy_precomputed():
    # Keep this specific random seed to match the precomputed Matlab result.
    rstate = np.random.RandomState(1234)
    # max RMSE (was 1.0e-10, is reduced to 5.0e-5 due to different coefficents)
    epsilon = 5.0e-5
    epsilon_pywt_coeffs = 1.0e-10
    for wavelet in wavelets:
        w = pywt.Wavelet(wavelet)
        for N in _get_data_sizes(w):
            data = rstate.randn(N)
            for pmode, mmode in modes:
                ma, md = _load_matlab_result(data, wavelet, mmode)
                yield _check_accuracy, data, w, pmode, ma, md, wavelet, epsilon
                ma, md = _load_matlab_result_pywt_coeffs(data, wavelet, mmode)
                yield _check_accuracy, data, w, pmode, ma, md, wavelet, epsilon_pywt_coeffs


def _compute_matlab_result(data, wavelet, mmode):
    """ Compute the result using MATLAB.

    This function assumes that the Matlab variables `wavelet` and `data` have
    already been set externally.
    """
    if np.any((wavelet == np.array(['coif6', 'coif7', 'coif8', 'coif9', 'coif10', 'coif11', 'coif12', 'coif13', 'coif14', 'coif15', 'coif16', 'coif17'])),axis=0):
        w = pywt.Wavelet(wavelet)
        mlab.set_variable('Lo_D', w.dec_lo)
        mlab.set_variable('Hi_D', w.dec_hi)
        mlab_code = ("[ma, md] = dwt(data, Lo_D, Hi_D, 'mode', '%s');" % mmode)
    else:
        mlab_code = "[ma, md] = dwt(data, wavelet, 'mode', '%s');" % mmode
    res = mlab.run_code(mlab_code)
    if not res['success']:
        raise RuntimeError("Matlab failed to execute the provided code. "
                           "Check that the wavelet toolbox is installed.")
    # need np.asarray because sometimes the output is a single float64
    ma = np.asarray(mlab.get_variable('ma'))
    md = np.asarray(mlab.get_variable('md'))
    return ma, md


def _load_matlab_result(data, wavelet, mmode):
    """ Load the precomputed result.
    """
    N = len(data)
    ma_key = '_'.join([mmode, wavelet, str(N), 'ma'])
    md_key = '_'.join([mmode, wavelet, str(N), 'md'])
    if (ma_key not in matlab_result_dict) or \
            (md_key not in matlab_result_dict):
        raise KeyError(
            "Precompted Matlab result not found for wavelet: "
            "{0}, mode: {1}, size: {2}".format(wavelet, mmode, N))
    ma = matlab_result_dict[ma_key]
    md = matlab_result_dict[md_key]
    return ma, md


def _load_matlab_result_pywt_coeffs(data, wavelet, mmode):
    """ Load the precomputed result.
    """
    N = len(data)
    ma_key = '_'.join([mmode, wavelet, str(N), 'ma_pywtCoeffs'])
    md_key = '_'.join([mmode, wavelet, str(N), 'md_pywtCoeffs'])
    if (ma_key not in matlab_result_dict) or \
            (md_key not in matlab_result_dict):
        raise KeyError(
            "Precompted Matlab result not found for wavelet: "
            "{0}, mode: {1}, size: {2}".format(wavelet, mmode, N))
    ma = matlab_result_dict[ma_key]
    md = matlab_result_dict[md_key]
    return ma, md


def _check_accuracy(data, w, pmode, ma, md, wavelet, epsilon):
    # PyWavelets result
    pa, pd = pywt.dwt(data, w, pmode)

    # calculate error measures
    rms_a = np.sqrt(np.mean((pa - ma) ** 2))
    rms_d = np.sqrt(np.mean((pd - md) ** 2))

    msg = ('[RMS_A > EPSILON] for Mode: %s, Wavelet: %s, '
           'Length: %d, rms=%.3g' % (pmode, wavelet, len(data), rms_a))
    assert_(rms_a < epsilon, msg=msg)

    msg = ('[RMS_D > EPSILON] for Mode: %s, Wavelet: %s, '
           'Length: %d, rms=%.3g' % (pmode, wavelet, len(data), rms_d))
    assert_(rms_d < epsilon, msg=msg)


if __name__ == '__main__':
    run_module_suite()
