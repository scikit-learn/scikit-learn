# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
from scipy.fftpack import fft, ifft, fftfreq
from scipy.interpolate import interp1d
from ._warps_cy import _warp_fast
from ._radon_transform import sart_projection_update
from warnings import warn


__all__ = ['radon', 'order_angles_golden_ratio', 'iradon', 'iradon_sart']


def radon(image, theta=None, circle=None):
    """
    Calculates the radon transform of an image given specified
    projection angles.

    Parameters
    ----------
    image : array_like, dtype=float
        Input image. The rotation axis will be located in the pixel with
        indices ``(image.shape[0] // 2, image.shape[1] // 2)``.
    theta : array_like, dtype=float, optional (default np.arange(180))
        Projection angles (in degrees).
    circle : boolean, optional
        Assume image is zero outside the inscribed circle, making the
        width of each projection (the first dimension of the sinogram)
        equal to ``min(image.shape)``.
        The default behavior (None) is equivalent to False.

    Returns
    -------
    radon_image : ndarray
        Radon transform (sinogram).  The tomography rotation axis will lie
        at the pixel index ``radon_image.shape[0] // 2`` along the 0th
        dimension of ``radon_image``.

    References
    ----------
    .. [1] AC Kak, M Slaney, "Principles of Computerized Tomographic
           Imaging", IEEE Press 1988.
    .. [2] B.R. Ramesh, N. Srinivasa, K. Rajgopal, "An Algorithm for Computing
           the Discrete Radon Transform With Some Applications", Proceedings of
           the Fourth IEEE Region 10 International Conference, TENCON '89, 1989

    Notes
    -----
    Based on code of Justin K. Romberg
    (http://www.clear.rice.edu/elec431/projects96/DSP/bpanalysis.html)

    """
    if image.ndim != 2:
        raise ValueError('The input image must be 2-D')
    if theta is None:
        theta = np.arange(180)
    if circle is None:
        warn('The default of `circle` in `skimage.transform.radon` '
             'will change to `True` in version 0.15.')
        circle = False

    if circle:
        radius = min(image.shape) // 2
        c0, c1 = np.ogrid[0:image.shape[0], 0:image.shape[1]]
        reconstruction_circle = ((c0 - image.shape[0] // 2) ** 2
                                 + (c1 - image.shape[1] // 2) ** 2)
        reconstruction_circle = reconstruction_circle <= radius ** 2
        if not np.all(reconstruction_circle | (image == 0)):
            warn('Radon transform: image must be zero outside the '
                 'reconstruction circle')
        # Crop image to make it square
        slices = []
        for d in (0, 1):
            if image.shape[d] > min(image.shape):
                excess = image.shape[d] - min(image.shape)
                slices.append(slice(int(np.ceil(excess / 2)),
                                    int(np.ceil(excess / 2)
                                        + min(image.shape))))
            else:
                slices.append(slice(None))
        slices = tuple(slices)
        padded_image = image[slices]
    else:
        diagonal = np.sqrt(2) * max(image.shape)
        pad = [int(np.ceil(diagonal - s)) for s in image.shape]
        new_center = [(s + p) // 2 for s, p in zip(image.shape, pad)]
        old_center = [s // 2 for s in image.shape]
        pad_before = [nc - oc for oc, nc in zip(old_center, new_center)]
        pad_width = [(pb, p - pb) for pb, p in zip(pad_before, pad)]
        padded_image = np.pad(image, pad_width, mode='constant',
                              constant_values=0)
    # padded_image is always square
    assert padded_image.shape[0] == padded_image.shape[1]
    radon_image = np.zeros((padded_image.shape[0], len(theta)))
    center = padded_image.shape[0] // 2

    shift0 = np.array([[1, 0, -center],
                       [0, 1, -center],
                       [0, 0, 1]])
    shift1 = np.array([[1, 0, center],
                       [0, 1, center],
                       [0, 0, 1]])

    def build_rotation(theta):
        T = np.deg2rad(theta)
        R = np.array([[np.cos(T), np.sin(T), 0],
                      [-np.sin(T), np.cos(T), 0],
                      [0, 0, 1]])
        return shift1.dot(R).dot(shift0)

    for i in range(len(theta)):
        rotated = _warp_fast(padded_image, build_rotation(theta[i]))
        radon_image[:, i] = rotated.sum(0)
    return radon_image


def _sinogram_circle_to_square(sinogram):
    diagonal = int(np.ceil(np.sqrt(2) * sinogram.shape[0]))
    pad = diagonal - sinogram.shape[0]
    old_center = sinogram.shape[0] // 2
    new_center = diagonal // 2
    pad_before = new_center - old_center
    pad_width = ((pad_before, pad - pad_before), (0, 0))
    return np.pad(sinogram, pad_width, mode='constant', constant_values=0)


def iradon(radon_image, theta=None, output_size=None,
           filter="ramp", interpolation="linear", circle=None):
    """
    Inverse radon transform.

    Reconstruct an image from the radon transform, using the filtered
    back projection algorithm.

    Parameters
    ----------
    radon_image : array_like, dtype=float
        Image containing radon transform (sinogram). Each column of
        the image corresponds to a projection along a different angle. The
        tomography rotation axis should lie at the pixel index
        ``radon_image.shape[0] // 2`` along the 0th dimension of
        ``radon_image``.
    theta : array_like, dtype=float, optional
        Reconstruction angles (in degrees). Default: m angles evenly spaced
        between 0 and 180 (if the shape of `radon_image` is (N, M)).
    output_size : int
        Number of rows and columns in the reconstruction.
    filter : str, optional (default ramp)
        Filter used in frequency domain filtering. Ramp filter used by default.
        Filters available: ramp, shepp-logan, cosine, hamming, hann.
        Assign None to use no filter.
    interpolation : str, optional (default 'linear')
        Interpolation method used in reconstruction. Methods available:
        'linear', 'nearest', and 'cubic' ('cubic' is slow).
    circle : boolean, optional
        Assume the reconstructed image is zero outside the inscribed circle.
        Also changes the default output_size to match the behaviour of
        ``radon`` called with ``circle=True``.
        The default behavior (None) is equivalent to False.

    Returns
    -------
    reconstructed : ndarray
        Reconstructed image. The rotation axis will be located in the pixel
        with indices
        ``(reconstructed.shape[0] // 2, reconstructed.shape[1] // 2)``.

    References
    ----------
    .. [1] AC Kak, M Slaney, "Principles of Computerized Tomographic
           Imaging", IEEE Press 1988.
    .. [2] B.R. Ramesh, N. Srinivasa, K. Rajgopal, "An Algorithm for Computing
           the Discrete Radon Transform With Some Applications", Proceedings of
           the Fourth IEEE Region 10 International Conference, TENCON '89, 1989

    Notes
    -----
    It applies the Fourier slice theorem to reconstruct an image by
    multiplying the frequency domain of the filter with the FFT of the
    projection data. This algorithm is called filtered back projection.

    """
    if radon_image.ndim != 2:
        raise ValueError('The input image must be 2-D')
    if theta is None:
        m, n = radon_image.shape
        theta = np.linspace(0, 180, n, endpoint=False)
    else:
        theta = np.asarray(theta)
    if len(theta) != radon_image.shape[1]:
        raise ValueError("The given ``theta`` does not match the number of "
                         "projections in ``radon_image``.")
    interpolation_types = ('linear', 'nearest', 'cubic')
    if interpolation not in interpolation_types:
        raise ValueError("Unknown interpolation: %s" % interpolation)
    if not output_size:
        # If output size not specified, estimate from input radon image
        if circle:
            output_size = radon_image.shape[0]
        else:
            output_size = int(np.floor(np.sqrt((radon_image.shape[0]) ** 2
                                               / 2.0)))
    if circle is None:
        warn('The default of `circle` in `skimage.transform.iradon` '
             'will change to `True` in version 0.15.')
        circle = False
    if circle:
        radon_image = _sinogram_circle_to_square(radon_image)

    th = (np.pi / 180.0) * theta
    # resize image to next power of two (but no less than 64) for
    # Fourier analysis; speeds up Fourier and lessens artifacts
    projection_size_padded = \
        max(64, int(2 ** np.ceil(np.log2(2 * radon_image.shape[0]))))
    pad_width = ((0, projection_size_padded - radon_image.shape[0]), (0, 0))
    img = np.pad(radon_image, pad_width, mode='constant', constant_values=0)

    # Construct the Fourier filter
    f = fftfreq(projection_size_padded).reshape(-1, 1)   # digital frequency
    omega = 2 * np.pi * f                                # angular frequency
    fourier_filter = 2 * np.abs(f)                       # ramp filter
    if filter == "ramp":
        pass
    elif filter == "shepp-logan":
        # Start from first element to avoid divide by zero
        fourier_filter[1:] = fourier_filter[1:] * np.sin(omega[1:]) / omega[1:]
    elif filter == "cosine":
        fourier_filter *= np.cos(omega)
    elif filter == "hamming":
        fourier_filter *= (0.54 + 0.46 * np.cos(omega / 2))
    elif filter == "hann":
        fourier_filter *= (1 + np.cos(omega / 2)) / 2
    elif filter is None:
        fourier_filter[:] = 1
    else:
        raise ValueError("Unknown filter: %s" % filter)
    # Apply filter in Fourier domain
    projection = fft(img, axis=0) * fourier_filter
    radon_filtered = np.real(ifft(projection, axis=0))

    # Resize filtered image back to original size
    radon_filtered = radon_filtered[:radon_image.shape[0], :]
    reconstructed = np.zeros((output_size, output_size))
    # Determine the center of the projections (= center of sinogram)
    mid_index = radon_image.shape[0] // 2

    [X, Y] = np.mgrid[0:output_size, 0:output_size]
    xpr = X - int(output_size) // 2
    ypr = Y - int(output_size) // 2

    # Reconstruct image by interpolation
    for i in range(len(theta)):
        t = ypr * np.cos(th[i]) - xpr * np.sin(th[i])
        x = np.arange(radon_filtered.shape[0]) - mid_index
        if interpolation == 'linear':
            backprojected = np.interp(t, x, radon_filtered[:, i],
                                      left=0, right=0)
        else:
            interpolant = interp1d(x, radon_filtered[:, i], kind=interpolation,
                                   bounds_error=False, fill_value=0)
            backprojected = interpolant(t)
        reconstructed += backprojected
    if circle:
        radius = output_size // 2
        reconstruction_circle = (xpr ** 2 + ypr ** 2) <= radius ** 2
        reconstructed[~reconstruction_circle] = 0.

    return reconstructed * np.pi / (2 * len(th))


def order_angles_golden_ratio(theta):
    """
    Order angles to reduce the amount of correlated information
    in subsequent projections.

    Parameters
    ----------
    theta : 1D array of floats
        Projection angles in degrees. Duplicate angles are not allowed.

    Returns
    -------
    indices_generator : generator yielding unsigned integers
        The returned generator yields indices into ``theta`` such that
        ``theta[indices]`` gives the approximate golden ratio ordering
        of the projections. In total, ``len(theta)`` indices are yielded.
        All non-negative integers < ``len(theta)`` are yielded exactly once.

    Notes
    -----
    The method used here is that of the golden ratio introduced
    by T. Kohler.

    References
    ----------
    .. [1] Kohler, T. "A projection access scheme for iterative
           reconstruction based on the golden section." Nuclear Science
           Symposium Conference Record, 2004 IEEE. Vol. 6. IEEE, 2004.
    .. [2] Winkelmann, Stefanie, et al. "An optimal radial profile order
           based on the Golden Ratio for time-resolved MRI."
           Medical Imaging, IEEE Transactions on 26.1 (2007): 68-76.
    """
    interval = 180

    def angle_distance(a, b):
        difference = a - b
        return min(abs(difference % interval), abs(difference % -interval))

    remaining = list(np.argsort(theta))   # indices into theta
    # yield an arbitrary angle to start things off
    index = remaining.pop(0)
    angle = theta[index]
    yield index
    # determine subsequent angles using the golden ratio method
    angle_increment = interval * (1 - (np.sqrt(5) - 1) / 2)
    while remaining:
        angle = (angle + angle_increment) % interval
        insert_point = np.searchsorted(theta[remaining], angle)
        index_below = insert_point - 1
        index_above = 0 if insert_point == len(remaining) else insert_point
        distance_below = angle_distance(angle, theta[remaining[index_below]])
        distance_above = angle_distance(angle, theta[remaining[index_above]])
        if distance_below < distance_above:
            yield remaining.pop(index_below)
        else:
            yield remaining.pop(index_above)


def iradon_sart(radon_image, theta=None, image=None, projection_shifts=None,
                clip=None, relaxation=0.15):
    """
    Inverse radon transform

    Reconstruct an image from the radon transform, using a single iteration of
    the Simultaneous Algebraic Reconstruction Technique (SART) algorithm.

    Parameters
    ----------
    radon_image : 2D array, dtype=float
        Image containing radon transform (sinogram). Each column of
        the image corresponds to a projection along a different angle. The
        tomography rotation axis should lie at the pixel index
        ``radon_image.shape[0] // 2`` along the 0th dimension of
        ``radon_image``.
    theta : 1D array, dtype=float, optional
        Reconstruction angles (in degrees). Default: m angles evenly spaced
        between 0 and 180 (if the shape of `radon_image` is (N, M)).
    image : 2D array, dtype=float, optional
        Image containing an initial reconstruction estimate. Shape of this
        array should be ``(radon_image.shape[0], radon_image.shape[0])``. The
        default is an array of zeros.
    projection_shifts : 1D array, dtype=float
        Shift the projections contained in ``radon_image`` (the sinogram) by
        this many pixels before reconstructing the image. The i'th value
        defines the shift of the i'th column of ``radon_image``.
    clip : length-2 sequence of floats
        Force all values in the reconstructed tomogram to lie in the range
        ``[clip[0], clip[1]]``
    relaxation : float
        Relaxation parameter for the update step. A higher value can
        improve the convergence rate, but one runs the risk of instabilities.
        Values close to or higher than 1 are not recommended.

    Returns
    -------
    reconstructed : ndarray
        Reconstructed image. The rotation axis will be located in the pixel
        with indices
        ``(reconstructed.shape[0] // 2, reconstructed.shape[1] // 2)``.

    Notes
    -----
    Algebraic Reconstruction Techniques are based on formulating the tomography
    reconstruction problem as a set of linear equations. Along each ray,
    the projected value is the sum of all the values of the cross section along
    the ray. A typical feature of SART (and a few other variants of algebraic
    techniques) is that it samples the cross section at equidistant points
    along the ray, using linear interpolation between the pixel values of the
    cross section. The resulting set of linear equations are then solved using
    a slightly modified Kaczmarz method.

    When using SART, a single iteration is usually sufficient to obtain a good
    reconstruction. Further iterations will tend to enhance high-frequency
    information, but will also often increase the noise.

    References
    ----------
    .. [1] AC Kak, M Slaney, "Principles of Computerized Tomographic
           Imaging", IEEE Press 1988.
    .. [2] AH Andersen, AC Kak, "Simultaneous algebraic reconstruction
           technique (SART): a superior implementation of the ART algorithm",
           Ultrasonic Imaging 6 pp 81--94 (1984)
    .. [3] S Kaczmarz, "Angenäherte auflösung von systemen linearer
           gleichungen", Bulletin International de l’Academie Polonaise des
           Sciences et des Lettres 35 pp 355--357 (1937)
    .. [4] Kohler, T. "A projection access scheme for iterative
           reconstruction based on the golden section." Nuclear Science
           Symposium Conference Record, 2004 IEEE. Vol. 6. IEEE, 2004.
    .. [5] Kaczmarz' method, Wikipedia,
           http://en.wikipedia.org/wiki/Kaczmarz_method
    """
    if radon_image.ndim != 2:
        raise ValueError('radon_image must be two dimensional')
    reconstructed_shape = (radon_image.shape[0], radon_image.shape[0])
    if theta is None:
        theta = np.linspace(0, 180, radon_image.shape[1], endpoint=False)
    elif theta.shape != (radon_image.shape[1],):
        raise ValueError('Shape of theta (%s) does not match the '
                         'number of projections (%d)'
                         % (projection_shifts.shape, radon_image.shape[1]))
    if image is None:
        image = np.zeros(reconstructed_shape, dtype=np.float)
    elif image.shape != reconstructed_shape:
        raise ValueError('Shape of image (%s) does not match first dimension '
                         'of radon_image (%s)'
                         % (image.shape, reconstructed_shape))
    if projection_shifts is None:
        projection_shifts = np.zeros((radon_image.shape[1],), dtype=np.float)
    elif projection_shifts.shape != (radon_image.shape[1],):
        raise ValueError('Shape of projection_shifts (%s) does not match the '
                         'number of projections (%d)'
                         % (projection_shifts.shape, radon_image.shape[1]))
    if clip is not None:
        if len(clip) != 2:
            raise ValueError('clip must be a length-2 sequence')
        clip = (float(clip[0]), float(clip[1]))
    relaxation = float(relaxation)

    for angle_index in order_angles_golden_ratio(theta):
        image_update = sart_projection_update(image, theta[angle_index],
                                              radon_image[:, angle_index],
                                              projection_shifts[angle_index])
        image += relaxation * image_update
        if clip is not None:
            image = np.clip(image, clip[0], clip[1])
    return image
