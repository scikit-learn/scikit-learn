"""
Benchmarks of Non-Negative Matrix Factorization
"""
# Authors: Tom Dupre la Tour (benchmark)
#          Chih-Jen Linn (original projected gradient NMF implementation)
#          Anthony Di Franco (projected gradient, Python and NumPy port)
# License: BSD 3 clause

from __future__ import print_function
from time import time
import sys
import warnings
import numbers

import numpy as np
import matplotlib.pyplot as plt
import pandas

from sklearn.utils.testing import ignore_warnings
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.decomposition.nmf import NMF
from sklearn.decomposition.nmf import _initialize_nmf
from sklearn.decomposition.nmf import _beta_divergence
from sklearn.decomposition.nmf import INTEGER_TYPES, _check_init
from sklearn.utils import Memory
from sklearn.exceptions import ConvergenceWarning
from sklearn.utils.extmath import safe_sparse_dot, squared_norm
from sklearn.utils import check_array
from sklearn.utils.validation import check_is_fitted, check_non_negative


mem = Memory(cachedir='.', verbose=0)

###################
# Start of _PGNMF #
###################
# This class implements a projected gradient solver for the NMF.
# The projected gradient solver was removed from scikit-learn in version 0.19,
# and a simplified copy is used here for comparison purpose only.
# It is not tested, and it may change or disappear without notice.


def _norm(x):
    """Dot product-based Euclidean norm implementation
    See: http://fseoane.net/blog/2011/computing-the-vector-norm/
    """
    return np.sqrt(squared_norm(x))


def _nls_subproblem(X, W, H, tol, max_iter, alpha=0., l1_ratio=0.,
                    sigma=0.01, beta=0.1):
    """Non-negative least square solver
    Solves a non-negative least squares subproblem using the projected
    gradient descent algorithm.
    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
        Constant matrix.
    W : array-like, shape (n_samples, n_components)
        Constant matrix.
    H : array-like, shape (n_components, n_features)
        Initial guess for the solution.
    tol : float
        Tolerance of the stopping condition.
    max_iter : int
        Maximum number of iterations before timing out.
    alpha : double, default: 0.
        Constant that multiplies the regularization terms. Set it to zero to
        have no regularization.
    l1_ratio : double, default: 0.
        The regularization mixing parameter, with 0 <= l1_ratio <= 1.
        For l1_ratio = 0 the penalty is an L2 penalty.
        For l1_ratio = 1 it is an L1 penalty.
        For 0 < l1_ratio < 1, the penalty is a combination of L1 and L2.
    sigma : float
        Constant used in the sufficient decrease condition checked by the line
        search.  Smaller values lead to a looser sufficient decrease condition,
        thus reducing the time taken by the line search, but potentially
        increasing the number of iterations of the projected gradient
        procedure. 0.01 is a commonly used value in the optimization
        literature.
    beta : float
        Factor by which the step size is decreased (resp. increased) until
        (resp. as long as) the sufficient decrease condition is satisfied.
        Larger values allow to find a better step size but lead to longer line
        search. 0.1 is a commonly used value in the optimization literature.
    Returns
    -------
    H : array-like, shape (n_components, n_features)
        Solution to the non-negative least squares problem.
    grad : array-like, shape (n_components, n_features)
        The gradient.
    n_iter : int
        The number of iterations done by the algorithm.
    References
    ----------
    C.-J. Lin. Projected gradient methods for non-negative matrix
    factorization. Neural Computation, 19(2007), 2756-2779.
    http://www.csie.ntu.edu.tw/~cjlin/nmf/
    """
    WtX = safe_sparse_dot(W.T, X)
    WtW = np.dot(W.T, W)

    # values justified in the paper (alpha is renamed gamma)
    gamma = 1
    for n_iter in range(1, max_iter + 1):
        grad = np.dot(WtW, H) - WtX
        if alpha > 0 and l1_ratio == 1.:
            grad += alpha
        elif alpha > 0:
            grad += alpha * (l1_ratio + (1 - l1_ratio) * H)

        # The following multiplication with a boolean array is more than twice
        # as fast as indexing into grad.
        if _norm(grad * np.logical_or(grad < 0, H > 0)) < tol:
            break

        Hp = H

        for inner_iter in range(20):
            # Gradient step.
            Hn = H - gamma * grad
            # Projection step.
            Hn *= Hn > 0
            d = Hn - H
            gradd = np.dot(grad.ravel(), d.ravel())
            dQd = np.dot(np.dot(WtW, d).ravel(), d.ravel())
            suff_decr = (1 - sigma) * gradd + 0.5 * dQd < 0
            if inner_iter == 0:
                decr_gamma = not suff_decr

            if decr_gamma:
                if suff_decr:
                    H = Hn
                    break
                else:
                    gamma *= beta
            elif not suff_decr or (Hp == Hn).all():
                H = Hp
                break
            else:
                gamma /= beta
                Hp = Hn

    if n_iter == max_iter:
        warnings.warn("Iteration limit reached in nls subproblem.",
                      ConvergenceWarning)

    return H, grad, n_iter


def _fit_projected_gradient(X, W, H, tol, max_iter, nls_max_iter, alpha,
                            l1_ratio):
    gradW = (np.dot(W, np.dot(H, H.T)) -
             safe_sparse_dot(X, H.T, dense_output=True))
    gradH = (np.dot(np.dot(W.T, W), H) -
             safe_sparse_dot(W.T, X, dense_output=True))

    init_grad = squared_norm(gradW) + squared_norm(gradH.T)
    # max(0.001, tol) to force alternating minimizations of W and H
    tolW = max(0.001, tol) * np.sqrt(init_grad)
    tolH = tolW

    for n_iter in range(1, max_iter + 1):
        # stopping condition as discussed in paper
        proj_grad_W = squared_norm(gradW * np.logical_or(gradW < 0, W > 0))
        proj_grad_H = squared_norm(gradH * np.logical_or(gradH < 0, H > 0))

        if (proj_grad_W + proj_grad_H) / init_grad < tol ** 2:
            break

        # update W
        Wt, gradWt, iterW = _nls_subproblem(X.T, H.T, W.T, tolW, nls_max_iter,
                                            alpha=alpha, l1_ratio=l1_ratio)
        W, gradW = Wt.T, gradWt.T

        if iterW == 1:
            tolW = 0.1 * tolW

        # update H
        H, gradH, iterH = _nls_subproblem(X, W, H, tolH, nls_max_iter,
                                          alpha=alpha, l1_ratio=l1_ratio)
        if iterH == 1:
            tolH = 0.1 * tolH

    H[H == 0] = 0   # fix up negative zeros

    if n_iter == max_iter:
        Wt, _, _ = _nls_subproblem(X.T, H.T, W.T, tolW, nls_max_iter,
                                   alpha=alpha, l1_ratio=l1_ratio)
        W = Wt.T

    return W, H, n_iter


class _PGNMF(NMF):
    """Non-Negative Matrix Factorization (NMF) with projected gradient solver.

    This class is private and for comparison purpose only.
    It may change or disappear without notice.

    """
    def __init__(self, n_components=None, solver='pg', init=None,
                 tol=1e-4, max_iter=200, random_state=None,
                 alpha=0., l1_ratio=0., nls_max_iter=10):
        super(_PGNMF, self).__init__(
            n_components=n_components, init=init, solver=solver, tol=tol,
            max_iter=max_iter, random_state=random_state, alpha=alpha,
            l1_ratio=l1_ratio)
        self.nls_max_iter = nls_max_iter

    def fit(self, X, y=None, **params):
        self.fit_transform(X, **params)
        return self

    def transform(self, X):
        check_is_fitted(self, 'components_')
        H = self.components_
        W, _, self.n_iter_ = self._fit_transform(X, H=H, update_H=False)
        return W

    def inverse_transform(self, W):
        check_is_fitted(self, 'components_')
        return np.dot(W, self.components_)

    def fit_transform(self, X, y=None, W=None, H=None):
        W, H, self.n_iter = self._fit_transform(X, W=W, H=H, update_H=True)
        self.components_ = H
        return W

    def _fit_transform(self, X, y=None, W=None, H=None, update_H=True):
        X = check_array(X, accept_sparse=('csr', 'csc'))
        check_non_negative(X, "NMF (input X)")

        n_samples, n_features = X.shape
        n_components = self.n_components
        if n_components is None:
            n_components = n_features

        if (not isinstance(n_components, INTEGER_TYPES) or
                n_components <= 0):
            raise ValueError("Number of components must be a positive integer;"
                             " got (n_components=%r)" % n_components)
        if not isinstance(self.max_iter, INTEGER_TYPES) or self.max_iter < 0:
            raise ValueError("Maximum number of iterations must be a positive "
                             "integer; got (max_iter=%r)" % self.max_iter)
        if not isinstance(self.tol, numbers.Number) or self.tol < 0:
            raise ValueError("Tolerance for stopping criteria must be "
                             "positive; got (tol=%r)" % self.tol)

        # check W and H, or initialize them
        if self.init == 'custom' and update_H:
            _check_init(H, (n_components, n_features), "NMF (input H)")
            _check_init(W, (n_samples, n_components), "NMF (input W)")
        elif not update_H:
            _check_init(H, (n_components, n_features), "NMF (input H)")
            W = np.zeros((n_samples, n_components))
        else:
            W, H = _initialize_nmf(X, n_components, init=self.init,
                                   random_state=self.random_state)

        if update_H:  # fit_transform
            W, H, n_iter = _fit_projected_gradient(
                X, W, H, self.tol, self.max_iter, self.nls_max_iter,
                self.alpha, self.l1_ratio)
        else:  # transform
            Wt, _, n_iter = _nls_subproblem(X.T, H.T, W.T, self.tol,
                                            self.nls_max_iter,
                                            alpha=self.alpha,
                                            l1_ratio=self.l1_ratio)
            W = Wt.T

        if n_iter == self.max_iter and self.tol > 0:
            warnings.warn("Maximum number of iteration %d reached. Increase it"
                          " to improve convergence." % self.max_iter,
                          ConvergenceWarning)

        return W, H, n_iter

#################
# End of _PGNMF #
#################


def plot_results(results_df, plot_name):
    if results_df is None:
        return None

    plt.figure(figsize=(16, 6))
    colors = 'bgr'
    markers = 'ovs'
    ax = plt.subplot(1, 3, 1)
    for i, init in enumerate(np.unique(results_df['init'])):
        plt.subplot(1, 3, i + 1, sharex=ax, sharey=ax)
        for j, method in enumerate(np.unique(results_df['method'])):
            mask = np.logical_and(results_df['init'] == init,
                                  results_df['method'] == method)
            selected_items = results_df[mask]

            plt.plot(selected_items['time'], selected_items['loss'],
                     color=colors[j % len(colors)], ls='-',
                     marker=markers[j % len(markers)],
                     label=method)

        plt.legend(loc=0, fontsize='x-small')
        plt.xlabel("Time (s)")
        plt.ylabel("loss")
        plt.title("%s" % init)
    plt.suptitle(plot_name, fontsize=16)


@ignore_warnings(category=ConvergenceWarning)
# use joblib to cache the results.
# X_shape is specified in arguments for avoiding hashing X
@mem.cache(ignore=['X', 'W0', 'H0'])
def bench_one(name, X, W0, H0, X_shape, clf_type, clf_params, init,
              n_components, random_state):
    W = W0.copy()
    H = H0.copy()

    clf = clf_type(**clf_params)
    st = time()
    W = clf.fit_transform(X, W=W, H=H)
    end = time()
    H = clf.components_

    this_loss = _beta_divergence(X, W, H, 2.0, True)
    duration = end - st
    return this_loss, duration


def run_bench(X, clfs, plot_name, n_components, tol, alpha, l1_ratio):
    start = time()
    results = []
    for name, clf_type, iter_range, clf_params in clfs:
        print("Training %s:" % name)
        for rs, init in enumerate(('nndsvd', 'nndsvdar', 'random')):
            print("    %s %s: " % (init, " " * (8 - len(init))), end="")
            W, H = _initialize_nmf(X, n_components, init, 1e-6, rs)

            for max_iter in iter_range:
                clf_params['alpha'] = alpha
                clf_params['l1_ratio'] = l1_ratio
                clf_params['max_iter'] = max_iter
                clf_params['tol'] = tol
                clf_params['random_state'] = rs
                clf_params['init'] = 'custom'
                clf_params['n_components'] = n_components

                this_loss, duration = bench_one(name, X, W, H, X.shape,
                                                clf_type, clf_params,
                                                init, n_components, rs)

                init_name = "init='%s'" % init
                results.append((name, this_loss, duration, init_name))
                # print("loss: %.6f, time: %.3f sec" % (this_loss, duration))
                print(".", end="")
                sys.stdout.flush()
            print(" ")

    # Use a panda dataframe to organize the results
    results_df = pandas.DataFrame(results,
                                  columns="method loss time init".split())
    print("Total time = %0.3f sec\n" % (time() - start))

    # plot the results
    plot_results(results_df, plot_name)
    return results_df


def load_20news():
    print("Loading 20 newsgroups dataset")
    print("-----------------------------")
    from sklearn.datasets import fetch_20newsgroups
    dataset = fetch_20newsgroups(shuffle=True, random_state=1,
                                 remove=('headers', 'footers', 'quotes'))
    vectorizer = TfidfVectorizer(max_df=0.95, min_df=2, stop_words='english')
    tfidf = vectorizer.fit_transform(dataset.data)
    return tfidf


def load_faces():
    print("Loading Olivetti face dataset")
    print("-----------------------------")
    from sklearn.datasets import fetch_olivetti_faces
    faces = fetch_olivetti_faces(shuffle=True)
    return faces.data


def build_clfs(cd_iters, pg_iters, mu_iters):
    clfs = [("Coordinate Descent", NMF, cd_iters, {'solver': 'cd'}),
            ("Projected Gradient", _PGNMF, pg_iters, {'solver': 'pg'}),
            ("Multiplicative Update", NMF, mu_iters, {'solver': 'mu'}),
            ]
    return clfs


if __name__ == '__main__':
    alpha = 0.
    l1_ratio = 0.5
    n_components = 10
    tol = 1e-15

    # first benchmark on 20 newsgroup dataset: sparse, shape(11314, 39116)
    plot_name = "20 Newsgroups sparse dataset"
    cd_iters = np.arange(1, 30)
    pg_iters = np.arange(1, 6)
    mu_iters = np.arange(1, 30)
    clfs = build_clfs(cd_iters, pg_iters, mu_iters)
    X_20news = load_20news()
    run_bench(X_20news, clfs, plot_name, n_components, tol, alpha, l1_ratio)

    # second benchmark on Olivetti faces dataset: dense, shape(400, 4096)
    plot_name = "Olivetti Faces dense dataset"
    cd_iters = np.arange(1, 30)
    pg_iters = np.arange(1, 12)
    mu_iters = np.arange(1, 30)
    clfs = build_clfs(cd_iters, pg_iters, mu_iters)
    X_faces = load_faces()
    run_bench(X_faces, clfs, plot_name, n_components, tol, alpha, l1_ratio,)

    plt.show()
