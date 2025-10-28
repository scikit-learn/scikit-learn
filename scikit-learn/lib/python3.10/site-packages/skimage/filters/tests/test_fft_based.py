import math

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal
import scipy.fft as fftmodule

from skimage._shared.utils import _supported_float_type
from skimage.data import astronaut, coins
from skimage.filters import butterworth
from skimage.filters._fft_based import _get_nd_butterworth_filter


def _fft_centered(x):
    return fftmodule.fftshift(fftmodule.fftn(fftmodule.fftshift(x)))


@pytest.mark.parametrize(
    'dtype', [np.float16, np.float32, np.float64, np.uint8, np.int32]
)
@pytest.mark.parametrize('squared_butterworth', [False, True])
def test_butterworth_2D_zeros_dtypes(dtype, squared_butterworth):
    im = np.zeros((4, 4), dtype=dtype)
    filtered = butterworth(im, squared_butterworth=squared_butterworth)
    assert filtered.shape == im.shape
    assert filtered.dtype == _supported_float_type(dtype)
    assert_array_equal(im, filtered)


@pytest.mark.parametrize('squared_butterworth', [False, True])
@pytest.mark.parametrize('high_pass', [False, True])
# order chosen large enough that lowpass stopband always approaches 0
@pytest.mark.parametrize('order', [6, 10])
@pytest.mark.parametrize('cutoff', [0.2, 0.3])
def test_butterworth_cutoff(cutoff, order, high_pass, squared_butterworth):
    wfilt = _get_nd_butterworth_filter(
        shape=(512, 512),
        factor=cutoff,
        order=order,
        high_pass=high_pass,
        real=False,
        squared_butterworth=squared_butterworth,
    )
    # select DC frequence on first axis to plot profile along a single axis
    wfilt_profile = np.abs(wfilt[0])

    # Empirical chosen to pass for order=6. Can use a smaller tolerance at
    # higher orders.
    tol = 0.3 / order

    # should have amplitude of ~1.0 in the center of the passband
    if high_pass:
        assert abs(wfilt_profile[wfilt_profile.size // 2] - 1.0) < tol
    else:
        assert abs(wfilt_profile[0] - 1.0) < tol

    # should be close to the expected amplitude at the cutoff frequency
    f_cutoff = int(cutoff * wfilt.shape[0])
    if squared_butterworth:
        # expect 0.5 at the cutoff
        assert abs(wfilt_profile[f_cutoff] - 0.5) < tol
    else:
        # expect 1/sqrt(2) at the cutoff
        assert abs(wfilt_profile[f_cutoff] - 1 / math.sqrt(2)) < tol


@pytest.mark.parametrize('cutoff', [-0.01, 0.51])
def test_butterworth_invalid_cutoff(cutoff):
    with pytest.raises(ValueError):
        butterworth(np.ones((4, 4)), cutoff_frequency_ratio=cutoff)


@pytest.mark.parametrize("high_pass", [True, False])
@pytest.mark.parametrize('squared_butterworth', [False, True])
def test_butterworth_2D(high_pass, squared_butterworth):
    # rough check of high-pass vs. low-pass behavior via relative energy

    # adjust specified order so that lowpass stopband approaches 0
    order = 3 if squared_butterworth else 6

    im = np.random.randn(64, 128)
    filtered = butterworth(
        im,
        cutoff_frequency_ratio=0.20,
        order=order,
        high_pass=high_pass,
        squared_butterworth=squared_butterworth,
    )

    # Compute the energy at the outer edges of the Fourier domain
    # before and after filtering.
    im_fft = _fft_centered(im)
    im_fft = np.real(im_fft * np.conj(im_fft))
    filtered_fft = _fft_centered(filtered)
    filtered_fft = np.real(filtered_fft * np.conj(filtered_fft))
    outer_mask = np.ones(im.shape, dtype=bool)
    outer_mask[4:-4, 4:-4] = 0
    abs_filt_outer = filtered_fft[outer_mask].mean()
    abs_im_outer = im_fft[outer_mask].mean()

    # Compute energy near the center of the Fourier domain
    inner_sl = tuple(slice(s // 2 - 4, s // 2 + 4) for s in im.shape)
    abs_filt_inner = filtered_fft[inner_sl].mean()
    abs_im_inner = im_fft[inner_sl].mean()

    if high_pass:
        assert abs_filt_outer > 0.9 * abs_im_outer
        assert abs_filt_inner < 0.1 * abs_im_inner
    else:
        assert abs_filt_outer < 0.1 * abs_im_outer
        assert abs_filt_inner > 0.9 * abs_im_inner


@pytest.mark.parametrize("high_pass", [True, False])
@pytest.mark.parametrize('dtype', [np.float32, np.float64])
@pytest.mark.parametrize('squared_butterworth', [False, True])
def test_butterworth_2D_realfft(high_pass, dtype, squared_butterworth):
    """Filtering a real-valued array is equivalent to filtering a
    complex-valued array where the imaginary part is zero.
    """
    im = np.random.randn(32, 64).astype(dtype)
    kwargs = dict(
        cutoff_frequency_ratio=0.20,
        high_pass=high_pass,
        squared_butterworth=squared_butterworth,
    )

    expected_dtype = _supported_float_type(im.dtype)
    filtered_real = butterworth(im, **kwargs)
    assert filtered_real.dtype == expected_dtype

    cplx_dtype = np.promote_types(im.dtype, np.complex64)
    filtered_cplx = butterworth(im.astype(cplx_dtype), **kwargs)
    assert filtered_cplx.real.dtype == expected_dtype

    if expected_dtype == np.float64:
        rtol = atol = 1e-13
    else:
        rtol = atol = 1e-5
    assert_allclose(filtered_real, filtered_cplx.real, rtol=rtol, atol=atol)


def test_butterworth_3D_zeros():
    im = np.zeros((3, 4, 5))
    filtered = butterworth(im)
    assert filtered.shape == im.shape
    assert_array_equal(im, filtered)


def test_butterworth_4D_zeros():
    im = np.zeros((3, 4, 5, 6))
    filtered = butterworth(im)
    assert filtered.shape == im.shape
    assert_array_equal(im, filtered)


@pytest.mark.parametrize(
    "chan, dtype",
    [(0, np.float64), (1, np.complex128), (2, np.uint8), (3, np.int64)],
)
def test_butterworth_4D_channel(chan, dtype):
    im = np.zeros((3, 4, 5, 6), dtype=dtype)
    filtered = butterworth(im, channel_axis=chan)
    assert filtered.shape == im.shape


def test_butterworth_correctness_bw():
    small = coins()[180:190, 260:270]
    filtered = butterworth(small, cutoff_frequency_ratio=0.2)
    # fmt: off
    correct = np.array(
        [
            [ 28.63019362, -17.69023786,  26.95346957,  20.57423019, -15.1933463 , -28.05828136, -35.25135674, -25.70376951, -43.37121955, -16.87688457],
            [  4.62077869,  36.5726672 ,  28.41926375, -22.86436829, -25.32375274, -19.94182623,  -2.9666164 ,   6.62250413,   3.55910886, -33.15358921],
            [ 25.00377084,  34.2948942 , -15.13862785, -15.34354183, -12.68722526,  12.82729905,   5.21622357,  11.41087761,  16.33690526, -50.39790969],
            [ 72.62331496, -14.7924709 , -22.14868895,  -7.47854864,   9.66784721,  24.37625693,  12.5479457 ,  -1.38194367,   2.40079497, -26.61141413],
            [ 21.85962078, -56.73932031, -14.82425429,   4.10524297, -19.16561768, -48.19021687,   5.0258744 ,  28.82432166,   0.66992097,   9.8378842 ],
            [-54.93417679,  -5.12037233,  19.2956981 ,  38.56431593,  27.95408908,  -3.53103389,  23.75329532,  -6.92615359,  -8.50823024,   7.05743093],
            [-30.51016624,  -9.99691211,  -7.1080672 ,  23.67643315,   1.61919726,  12.94103905, -29.08141699, -11.56937511,  22.70988847,  32.04063285],
            [ -7.51780937, -30.27899181,  -2.57782655,  -1.58947887,  -2.13564576, -11.34039302,   1.59165041,  14.39173421, -14.15148821,  -2.21664717],
            [ 14.81167298,  -3.75274782,  18.41459894,  15.80334075, -19.7477109 ,  -3.68619619,  -2.9513036 , -10.17325791,  18.32438702,  18.68003971],
            [-50.53430811,  12.14152989,  17.69341877,   9.1858496 ,  12.1470914 ,   1.45865179,  61.08961357,  29.76775029, -11.04603619,  24.18621404],
        ]
    )
    # fmt: on
    assert_allclose(filtered, correct)


def test_butterworth_correctness_rgb():
    small = astronaut()[135:145, 205:215]
    filtered = butterworth(
        small, cutoff_frequency_ratio=0.3, high_pass=True, channel_axis=-1
    )
    correct = np.array(
        [
            [
                [-5.30292781e-01, 2.17985072e00, 2.86622486e00],
                [6.39360740e00, 9.30643715e00, 8.61226660e00],
                [5.47978436e00, 1.03641402e01, 1.02940329e01],
                [8.88312002e00, 7.29681652e00, 8.16021235e00],
                [4.67693778e00, 6.33135145e-01, -2.51296407e00],
                [3.21039522e00, -9.14893931e-01, -2.11543661e00],
                [4.61985125e-02, -6.10682453e00, -1.72837650e00],
                [-4.59492989e00, -7.35887525e00, -1.03532871e01],
                [-3.67859542e00, -4.36371621e00, -5.67371459e00],
                [-4.38264080e00, -6.08362280e00, -9.20394882e00],
            ],
            [
                [-2.46191390e00, -2.10996960e00, -1.41287606e00],
                [4.87042304e-01, 4.70686760e-01, 2.90817746e00],
                [9.33095004e-01, -2.11867564e-01, 3.10917925e00],
                [-2.35660768e00, -1.35043153e00, -2.67062162e00],
                [-1.22363424e00, 1.11155488e-01, 1.25392954e00],
                [-1.05667680e00, 1.58195605e-01, 6.11873557e-01],
                [-4.12128910e00, -3.55994486e00, -8.75303054e00],
                [2.47171790e00, 2.70762582e00, 5.69543552e00],
                [6.97042504e-01, -2.24173305e00, 3.26477871e-01],
                [5.00195333e-01, -2.66024743e00, -1.87479563e00],
            ],
            [
                [-4.40136260e00, -4.02254309e00, -4.89246563e00],
                [-4.64563864e00, -6.21442755e00, -9.31399553e00],
                [-2.11532959e00, -2.58844609e00, -4.20629743e00],
                [-3.40862389e00, -3.29511853e00, -4.78220207e00],
                [-8.06768327e-01, -4.01651211e00, -2.84783939e00],
                [1.72379068e00, 1.00930709e00, 2.57575911e00],
                [-2.13771052e00, -1.75564076e00, -2.88676819e00],
                [-2.72015191e-01, -1.61610409e-01, 2.15906305e00],
                [-3.80356741e00, -7.30201675e-01, -3.79800352e00],
                [1.43534281e-01, -2.95121861e00, -2.67070135e00],
            ],
            [
                [1.03480271e00, 6.34545011e00, 3.53610283e00],
                [7.44740677e00, 9.97350707e00, 1.25152734e01],
                [3.10493189e00, 5.15553793e00, 6.48354940e00],
                [-7.89260096e-01, 3.04573015e-01, -1.43829810e00],
                [-1.46298411e00, 1.23095495e00, -1.33983509e00],
                [2.82986807e00, 2.80546223e00, 6.39492794e00],
                [-9.15293187e-01, 2.88688464e00, -9.69417480e-01],
                [4.50217964e00, 2.90410068e00, 5.39107589e00],
                [-5.71608069e-01, 1.78198962e00, -3.72062011e-01],
                [7.43870617e00, 8.78780364e00, 9.91142612e00],
            ],
            [
                [-1.01243616e01, -1.46725955e01, -1.63691866e01],
                [1.06512445e01, 7.84699418e00, 1.05428678e01],
                [4.75343829e00, 6.11329861e00, 2.81633365e00],
                [7.78936796e00, 9.63684277e00, 1.21495065e01],
                [5.19238043e00, 5.38165743e00, 8.03025884e00],
                [1.67424214e00, 2.25530135e00, 2.44161390e-01],
                [-3.18012002e-01, 1.99405335e00, -4.33960644e-01],
                [-1.21700957e00, -2.65973900e00, -6.31515766e-01],
                [-4.87805104e00, -5.55289609e00, -8.50052504e00],
                [1.43493808e01, 1.77252074e01, 1.92810954e01],
            ],
            [
                [-4.21712178e01, -4.47409535e01, -4.16375264e01],
                [1.18585381e00, 1.33732681e00, -1.45927283e00],
                [-4.83275742e00, -7.14344851e00, -7.59812923e00],
                [-7.13716513e00, -1.10025632e01, -1.16111397e01],
                [-5.00275373e00, -4.20410732e00, -4.93030043e00],
                [1.98421851e00, 2.68393141e00, 3.14898078e00],
                [1.97471502e00, -2.11937555e00, 2.04674150e00],
                [3.42916035e00, 4.98808524e00, 6.74436447e00],
                [-3.29035900e-01, -9.88239773e-01, 2.38909382e-02],
                [2.58646940e01, 2.40294324e01, 2.64535438e01],
            ],
            [
                [-2.72408341e01, -2.51637965e01, -2.95566971e01],
                [4.83667855e00, 5.63749968e00, 6.97366639e00],
                [6.18182884e00, 2.89519333e00, 6.15697112e00],
                [-1.03326540e00, -4.04702873e00, -5.50872246e00],
                [-6.92401355e00, -8.99374166e00, -9.43201766e00],
                [-7.24967366e00, -1.16225741e01, -1.26385982e01],
                [-7.79470220e00, -9.36143025e00, -7.13686778e00],
                [-1.22393834e01, -1.24094588e01, -1.70498337e01],
                [-4.65034860e00, -3.93071471e00, -4.65788605e00],
                [2.47461715e01, 2.46389560e01, 3.11843915e01],
            ],
            [
                [1.22510987e00, -3.30423292e00, -1.02785428e01],
                [-2.41285934e01, -1.60883211e01, -1.97064909e01],
                [9.77242420e00, 1.80847757e01, 2.01088714e01],
                [1.21335913e01, 1.16812821e01, 1.20531919e01],
                [-1.64052514e00, -4.45404672e00, -3.83103216e00],
                [9.97519545e-01, -6.32678881e00, -4.89843180e-01],
                [1.07464325e01, 2.62880708e01, 1.74691665e01],
                [2.31869805e00, 1.91935135e00, -9.65612833e-01],
                [-9.23026361e00, -1.75138706e01, -1.39243019e01],
                [4.60784836e00, 5.60845273e00, 5.28255564e00],
            ],
            [
                [9.68795318e00, 2.78851276e00, 5.45620353e00],
                [-1.02076906e01, -1.47926224e01, -1.43257049e01],
                [-2.17025353e00, 6.23446752e00, 5.21771748e00],
                [1.57074742e01, 2.17163634e01, 2.55600809e01],
                [7.03884578e00, 1.29273058e01, 7.50960315e00],
                [-6.69896692e00, -1.83433042e01, -1.60702492e01],
                [7.44877725e00, 1.28971365e01, 1.10234666e01],
                [5.25663607e00, 9.80648891e00, 1.22955858e01],
                [-7.44903684e00, -1.92670342e01, -1.68232131e01],
                [1.25609220e01, 2.09808909e01, 2.21425299e01],
            ],
            [
                [2.87825597e00, 3.37945227e00, 3.05777360e00],
                [1.18858884e00, -5.27430874e00, -6.96009863e00],
                [-7.55910235e00, -2.12068126e01, -2.06925790e01],
                [-1.47217788e01, -1.45626702e01, -1.56493571e01],
                [-5.60886203e00, 8.81908697e-01, 5.47367282e00],
                [-1.00478644e01, -8.01471176e00, -7.45670458e00],
                [3.61521638e00, 8.99194959e00, 4.93826323e00],
                [7.87025438e00, 1.34804191e01, 1.96899695e01],
                [-5.50012037e00, -6.40490471e00, -1.17265188e01],
                [6.17010624e00, 1.56199152e01, 1.79889524e01],
            ],
        ]
    )
    assert_allclose(filtered, correct)
