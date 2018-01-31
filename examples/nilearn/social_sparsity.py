"""
social sparsity: approximate overlapping group lasso

"""
# Author: Elvis Dohmatob, Gael Varoquaux
# License: simplified BSD

from math import sqrt, log
from numbers import Number
import numpy as np
from scipy import ndimage


def _fwhm2sigma(fwhm, voxel_size=None, affine=None):
    """Convert a FWHM value to sigma in a Gaussian kernel.
    """
    if voxel_size is None:
        if affine is None:
            raise ValueError(
                "voxel_size not provided, you must specify the affine")
        affine = affine[:3, :3]
        voxel_size = np.sqrt(np.sum(affine ** 2, axis=0))
    fwhm = np.asanyarray(fwhm)
    return fwhm / (sqrt(8. * log(2.)) * voxel_size)


def _prox_l21(img, alpha, grp_norms_squared=None, copy=False,
              tol=0.):
    # To avoid side effects, assign the raw img values on the side
    if copy:
        img = img.copy()
    if grp_norms_squared is None:
        grp_norms_squared = img ** 2
    grp_norms = np.sqrt(grp_norms_squared, out=grp_norms_squared)
    nz = (grp_norms > tol)
    img[nz] *= np.maximum(1. - alpha / grp_norms[nz], 0.)
    img[~nz] = 0.
    return img


def _gaussian_filternd(img, sigma, **kwargs):
    """Applies multi-dimensional Gaussian filter inplace."""
    if isinstance(sigma, Number):
        sigma = (sigma,) * img.ndim
    for axis, s in enumerate(sigma):
        ndimage.gaussian_filter1d(img, s, axis=axis, output=img, **kwargs)


def _grp_norms_squared(img, sigma, side_weights=.7, kernel="gaussian",
                       mode="constant", cval=0., n_resolutions=1):
    """Social sparsity as defined by eq 4 of Kowalski et al, 'Social Sparsity...'

    Parameters
    ----------
    side_weight: nonnegative float, optional (default 0.7)
        Weights of sides of neigborhood relative to center. A value of 1
        corresponds to the classical overlapping group-Lasso shrinkage
        operator.

    fwhm: int, optional (default 1)
        Size of neigbourhoods to consider, measured in mm's. This is a
        field-of-view (FOV) parameter and plays a rule similar to the fwhm in
        Gaussian kernels. The larger the radius, the smoother the prox.

    """
    if isinstance(sigma, Number):
        sigma = (sigma,) * img.ndim
    sigma = np.asanyarray(sigma)
    grp_norms_squared = img ** 2
    if kernel == "gaussian":
        _gaussian_filternd(grp_norms_squared, sigma, mode=mode, cval=cval)
    elif kernel == "pyramid":
        radius = np.floor(2 * sigma).astype(np.int)  # the "blur radius"
        diameter = 2 * radius + 1
        if side_weights == 1.:
            # use scipy's optimized rectangular uniform filter
            ndimage.uniform_filter(grp_norms_squared, size=diameter,
                                   output=grp_norms_squared, mode=mode,
                                   cval=cval)
        else:
            social_filter = np.full(diameter, 1.)
            social_filter *= side_weights

            # adjust weight at center of filter
            social_filter[np.ix_(*[[r] for r in radius])] = 1.
            # if img.ndim == 1:
            #     social_filter[radius] = 1.
            # elif img.ndim == 2:
            #     social_filter[radius[0], radius[1]] = 1.
            # elif img.ndim == 3:
            #     social_filter[radius[0], radius[1], radius[2]] = 1.
            # else:
            #     raise RuntimeError("WTF! img.ndim is %i." % img.ndim)

            # normalize filter weights to sum to 1
            social_filter /= social_filter.sum()

            # the actual convolution
            ndimage.filters.convolve(grp_norms_squared, social_filter,
                                     output=grp_norms_squared,
                                     mode="constant")
    else:
        raise ValueError("Unknown kernel: %s" % kernel)
    # else:
    #     # use ninja code from @gael
    #     grp_norms_squared = _neighboorhood_norms_squared(
    #         img, side_weights=side_weights)
    if n_resolutions > 1:
        grp_norms_squared += _grp_norms_squared(
            img, sigma, side_weights=side_weights, kernel=kernel,
            mode=mode, cval=cval, n_resolutions=n_resolutions - 1)
    return grp_norms_squared


def _prox_social_sparsity(img, alpha, fwhm, copy=False, affine=None,
                          voxel_size=None, side_weights=.7, kernel="gaussian",
                          mode="constant", cval=0., n_resolutions=1,
                          pos=False):
    """Social sparsity as defined by eq 4 of Kowalski et al, 'Social Sparsity...'

    Parameters
    ----------
    side_weight: nonnegative float, optional (default 0.7)
        Weights of sides of neigborhood relative to center. A value of 1
        corresponds to the classical overlapping group-Lasso shrinkage
        operator.

    fwhm: int, optional (default 1)
        Size of neigbourhoods to consider, measured in mm's. This is a
        field-of-view (FOV) parameter and plays a rule similar to the fwhm in
        Gaussian kernels. The larger the radius, the smoother the prox.

    """
    if copy:
        img = img.copy()
    sigma = _fwhm2sigma(fwhm, voxel_size=voxel_size, affine=affine)
    grp_norms_squared = _grp_norms_squared(
            img, sigma, side_weights=side_weights, kernel=kernel, mode=mode,
            cval=cval, n_resolutions=n_resolutions)
    img = _prox_l21(img, alpha, grp_norms_squared=grp_norms_squared,
                    copy=False)
    if pos:
        img = img.clip(min=0., out=img)
    return img


def _social_sparsity_alpha_grid(grad_loss, mask, fwhm, eps=1e-3,
                                affine=None, n_alphas=10, kernel="gaussian",
                                voxel_size=None, side_weights=.7,
                                mode="constant", cval=0., n_resolutions=1):
    """Computes grid of regularization parameters for social sparsity.

    Parameters
    ----------
    grad_loss: ndarray, shape (n_targets, n_features)
        Gradient of loss function at zero.

    mask: ndarray, shape (d1, d2, ...,)
        Contains n_features non-zero values.
    """
    imgs = np.zeros((len(grad_loss),) + mask.shape, dtype=grad_loss.dtype)
    imgs[:, mask] = grad_loss
    sigma = _fwhm2sigma(fwhm, voxel_size=voxel_size, affine=affine)
    grp_norms_squared  = [_grp_norms_squared(img, sigma,
                                             side_weights=side_weights,
                                             kernel=kernel, mode=mode,
                                             cval=cval,
                                             n_resolutions=n_resolutions)
                         for img in imgs]
    alpha_max = np.sqrt(np.max(grp_norms_squared))

    if n_alphas == 1:
        return np.array([alpha_max])
    alpha_min = alpha_max * eps
    return np.logspace(np.log10(alpha_min), np.log10(alpha_max),
                       num=n_alphas)[::-1]


if __name__ == "__main__":
    import matplotlib.pyplot as plt    
    from sklearn.utils import check_random_state    
    from nilearn.decoding.proximal_operators import _prox_tvl1
    from nilearn.plotting.cm import cold_hot

    rng = check_random_state(0)

    n_resolutions = 3
    fwhm = 4
    radius = _fwhm2sigma(fwhm, 1) * 2
    noise_sigma = .5
    mode = "constant"
    img_size = 50
    img_shape = (img_size,) * 2
    mask = np.ones(img_shape, dtype=np.bool)
    n_blobs = 100
    blob_size = 5
    img  = np.zeros(img_shape)
    for _ in range(n_blobs):
        blob_size = rng.binomial(blob_size, .8)
        sign = rng.choice([-1., 1.])
        amp = rng.choice(range(2, 11))
        weight = amp * sign
        i = rng.choice(img_shape[0])
        j = rng.choice(img_shape[1])
        img[i: i + blob_size, j:j + blob_size] = weight

    orig = img.copy()
    img += rng.randn(*img_shape) * noise_sigma

    proxes = {}
    orig2 = np.sum(orig ** 2)
    for kernel in ["tv", "gaussian", "pyramid"]:
        if kernel != "tv":
            name = "social(kernel=%s)" % kernel
            alphas = _social_sparsity_alpha_grid(
                img.ravel()[None, :], mask, fwhm, kernel=kernel, voxel_size=1,
                n_alphas=100, eps=1e-4, n_resolutions=n_resolutions)
            scores = []
            for alpha in alphas:
                prox = _prox_social_sparsity(img, alpha, fwhm, voxel_size=1,
                                             mode=mode, kernel=kernel,
                                             side_weights=.2, copy=True,
                                             n_resolutions=n_resolutions)
                score = np.sum((prox - orig) ** 2)
                scores.append(score)
            # alpha_best = alphas[np.argmin(scores)]
            # alpha = alpha_best
            alpha = alphas[0] * .4
            prox = _prox_social_sparsity(img, alpha, fwhm, voxel_size=1,
                                         mode=mode, kernel=kernel,
                                         side_weights=.2, copy=True,
                                         n_resolutions=n_resolutions)
            # prox = denoise_tv_chambolle(prox, weight=.2)
            plt.plot(np.log10(alphas / alphas[0]), scores)
            plt.title(name)
            plt.grid("on")
            # plt.axvline(np.log10(alpha_best / alphas[0]), linestyle="--",
            #             label="$\\alpha = \\alpha^*$")
            plt.xlabel("$\\log_{10}(\\alpha / \\alpha_0)$")
            plt.ylabel("SURE statistic")
            plt.show()
        else:
            name = kernel
            prox, _ = _prox_tvl1(img, l1_ratio=.5, weight=1, verbose=2)
        proxes[name] = prox
    ncols = 2 + len(proxes)
    _, axes = plt.subplots(1, ncols, figsize=(3 * ncols, 4))
    for which, ax, data in zip(["original", "noisy"] + list(proxes.keys()),
                               axes.ravel(),
                               [orig, img] + list(proxes.values())):
        im = ax.imshow(data, cmap=cold_hot)
        ax.axis("off")
        ax.set_title(which)
    # plt.colorbar(im)
    plt.tight_layout()
    plt.show()
