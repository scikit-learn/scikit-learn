from math import sqrt
import numpy as np
from scipy.optimize import check_grad
from sklearn.utils.testing import assert_less
from sklearn.utils import check_random_state
from sklearn.utils.testing import (assert_array_equal,
                                   assert_array_almost_equal)
from nilearn.decoding.objective_functions import _gradient, _div, _unmask
from nilearn.decoding.fista import mfista
from nilearn.decoding.proximal_operators import _prox_tvl1
from sklearn.linear_model.python_wrappers import _py_prox_l1, _py_proj_l1
from social_sparsity import _prox_social_sparsity
try:
    from modl.utils.math.enet import enet_projection
except ImportError:
    pass
from owl import prox_oscar, prox_owl
from ksupp import prox_ksupp


def _gram_schmidt(V, offset=0, normalize=True):
    """Computes modified Gram-Schmidt (MGS) orthogonalization of a set
    of vectors specified as rows of a 2D array V.

    Parameters
    ----------
    offset: int, optional (default 0)
    assumes that all vectors from index 0 to offset - 1 (inclusive) have
    already been processed.

    Returns
    -------
    V[offset:]
    """
    scales = {}
    for i in range(offset, len(V)):
        for j in range(i):
            if j not in scales:
                scales[j] = V[j].dot(V[j])
            if scales[j] > 0.:
                weight = V[j].dot(V[i]) / scales[j]
                V[i] -= weight * V[j]
        if normalize:
            scales[i] = V[i].dot(V[i])
            if scales[i] > 0.:
                V[i] /= sqrt(scales[i])
    return V[offset:]


def test_gram_schmidt():

    V = np.array([[3., 1.], [2., 2.]])
    assert_array_almost_equal(_gram_schmidt(V, normalize=False),
                              np.array([[3, 1], [-.4, 1.2]]))

    V = np.array([[1., 1., 1., 1.], [-1., 4., 4., 1.], [4., -2., 2., 0]])
    assert_array_almost_equal(_gram_schmidt(V, normalize=False),
                              np.array([[1., 1., 1., 1.], [-3., 2., 2., -1.],
                                        [1., -10. / 6., 7. / 3, -10. / 6.]]))

    V = np.array([[1., 1., 1., 1.], [-1., 4., 4., 1.], [4., -2., 2., 0]])
    V_ = V.copy()
    _gram_schmidt(V, offset=2, normalize=False)
    assert_array_equal(V[:2], V_[:2])


class ProximalOperator(object):
    """
    Solves for SSODL (Smooth Sparse Dictionary Learning) dictionary update.

    \argmin_{V} \frac{1}{2} .5/n * tr(VV^TA) - (1/n)tr(VB^T) +
                n * alpha \sum_j \varphi(v^j),

    where \varphi is a regularizer that imposes structure: sparsity and
    smoothness like GraphNet, TVL1, or social sparsity.

    References
    ==========
    [1] Dohmatob et al. "Learning brain regions via large-scale online
        structured sparse dictionary-learning", NIPS 2017
    [2] Varoquaux et al. "Social-sparsity brain decoders: faster spatial
        sparsity", PRNI 2016.
    [3] Kowalski et al. "Social sparsity! neighborhood systems enrich
        structured shrinkage operators", Transactions on Signal Processing
    """

    def __init__(self, which="enet", verbose=1, kernel="gaussian",
                 **params):
        self.which = which
        self.kernel = kernel
        self.params = params
        self.verbose = verbose

    def __call__(self, atom, weight, weight_scale=1., init=None):
        if weight_scale == 0.:
            atom[:] = 0.
            return atom
        n_voxels = len(atom)
        if weight_scale == 0.:
            atom[:] = 0.
            return atom
        pos = self.params.get("pos", True)
        l1_ratio = self.params.get("l1_ratio", 1.)
        variational_sparsity = self.params.get("variational_sparsity", False)
        if self.which == "n.o.p":
            return atom
        elif self.which == "enet":
            if l1_ratio != 1.:
                raise NotImplementedError("Non-lasso enet")
            if pos:
                atom[atom < 0.] = 0.
            out = np.zeros_like(atom)
            _py_proj_l1(atom, weight, weight_scale)
        elif self.which == "social":
            if weight_scale != 1.:
                atom /= weight_scale
                weight /= weight_scale
            params = {}
            for param in ["fwhm", "mask", "affine", "voxel_size",
                          "kernel"]:
                params[param] = self.params.get(param, None)
            side_weights = self.params.get("side_weights", .7)
            tmp = _prox_social_sparsity(
                _unmask(atom, params["mask"]), weight, params["fwhm"],
                affine=params["affine"], voxel_size=params["voxel_size"],
                pos=pos, kernel=self.kernel, side_weights=side_weights)[params["mask"]]
            atom[:] = tmp
        elif self.which == "gram-schmidt":
            raise NotImplementedError
            # dictionary[k] = _gram_schmidt(dictionary[:k + 1], offset=k)[-1]
        elif self.which == "enet variational":
            if l1_ratio != 1.:
                raise NotImplementedError("Non-lasso enet variational")
            _py_prox_l1(atom, weight, weight_scale)
        elif self.which == "ksupp":
            if weight_scale != 1.:
                atom /= weight_scale
                weight /= weight_scale
            out = prox_ksupp(atom, weight, verbose=0,
                             k=self.params.get("k", 3))
            atom = out[:]
        elif self.which == "owl":
            if weight_scale != 1.:
                atom /= weight_scale
                weight /= weight_scale
            sample_weights = self.params.get('sample_weights', 1.)
            atom[:] = prox_owl(atom, weight, sample_weights)
        elif self.which == "oscar":
            if weight_scale != 1.:
                atom /= weight_scale
                weight /= weight_scale
            alpha = self.params.get('alpha', .5 / len(atom))
            atom[:] = prox_oscar(atom, weight, alpha)
        elif self.which in ["tv-l1", "graph-net"]:
            if weight_scale != 1.:
                atom /= weight_scale
                weight = weight / weight_scale
            # misc
            max_iter = self.params.get("max_iter", 1000)
            tol = self.params.get("tol", 1e-3)
            mask = self.params["mask"]
            flat_mask = mask.ravel()
            l1_weight = weight * l1_ratio
            if self.which == "graph-net":
                spatial_weight = weight - l1_weight
            else:
                spatial_weight = 0.
            lap_lips = 4. * mask.ndim * spatial_weight
            loss_lips = 1.
            loss_lips *= 1.05
            lips = loss_lips + lap_lips

            def smooth_energy(v):
                """Smooth part of energy / cost function.
                """
                e = .5 * np.sum((v - atom) ** 2)
                if self.which == "graph-net":
                    lap = np.sum(_gradient(_unmask(v, mask))[:, mask] ** 2)
                    lap *= spatial_weight
                    lap *= .5
                    e += lap
                return e

            def nonsmooth_energy(v):
                """Non-smooth part of energy / cost function.
                """
                e = l1_weight * np.abs(v).sum()
                if self.which == "tv-l1":
                    gradient = _gradient(_unmask(v, mask))
                    tv_term = np.sum(np.sqrt(np.sum(gradient ** 2,
                                                    axis=0)))
                    tv_term *= spatial_weight
                    e += tv_term
                return e

            def total_energy(v):
                """Total energy / cost function.
                """
                return smooth_energy(v) + nonsmooth_energy(v)

            def grad(v):
                """Gradient of smooth part of energy / cost function.
                """
                grad = v - atom
                if self.which == "graph-net":
                    lap = -_div(_gradient(_unmask(v, mask)))[mask]
                    lap *= spatial_weight
                    grad += lap
                return grad

            def prox(v, stepsize, dgap_tol, init=None):
                """Proximal operator of non-smooth part of energy / cost function
                """
                if self.which == "graph-net":
                    info = dict(converged=True)
                    if variational_sparsity:
                        out = _prox_l1(v, stepsize * l1_weight, copy=False)
                    else:
                        out = np.ndarray(v.shape, dtype=v.dtype)
                        radius = self.params.get("radius", 1.)
                        enet_projection(v, out, radius, l1_ratio)
                    if pos:
                        out[out < 0.] = 0.
                elif self.which == "tv-l1":
                    v = _unmask(v, mask)
                    if init is not None:
                        init = _unmask(init, mask)
                    out, info = _prox_tvl1(v, init=init, l1_ratio=l1_ratio,
                                           weight=weight * stepsize,
                                           dgap_tol=dgap_tol, max_iter=100,
                                           verbose=self.verbose)
                    out = out.ravel()[flat_mask]
                else:
                    raise ValueError(
                        "Unknown value for self.which: %s" % (self.which))
                return out, info

            # for debugging
            check_gradient = self.params.get("check_grad", False)
            if check_gradient:
                rng = check_random_state(42)
                x0 = rng.randn(n_voxels)
                assert_less(check_grad(smooth_energy, grad, x0), 1e-3)

            # use FISTA update atom
            # init = dict(w=atom)
            check_lipschitz = self.params.get("check_lipschitz", False)
            out, _, _ = mfista(
                grad, prox, total_energy, lips, n_voxels, tol=tol,
                max_iter=max_iter, check_lipschitz=check_lipschitz,
                verbose=self.verbose,  # init=init
            )
            atom[:] = out
        else:
            raise NotImplementedError("which=%s" % self.which)

        return atom

    def __repr__(self):
        s = "%s(which='%s', weight=%s, verbose=%s, %s)" % (
            self.__class__.__name__, self.which, self.weight, self.verbose,
            ", ".join(["%s=%s" % (param, self.params[param])
                       for param in self.params]))
        return s
