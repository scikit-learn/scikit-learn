"""
Benchmarks of Lasso vs LassoLars

First, we fix a training set and increase the number of
samples. Then we plot the computation time as function of
the number of samples.

In the second benchmark, we increase the number of dimensions of the
training set. Then we plot the computation time as function of
the number of dimensions.

In both cases, only 10% of the features are informative.
"""
from time import time
import numpy as np
from sklearn.datasets.samples_generator import make_regression
from sklearn.utils import check_random_state
from sklearn.linear_model import cd_fast, cd_fast2


def compute_bench(alpha, n_samples, n_features, precompute,
                  n_targets=5, max_iter=1000, tol=1e-4, l1_ratio=0.,
                  backend="legacy", multi_task=False, random_state=0):
    it = 0
    results = []
    rng = check_random_state(random_state)
    for ns in n_samples:
        for nf in n_features:
            it += 1
            print('==================')
            print('Iteration %s of %s' % (it, max(len(n_samples),
                                          len(n_features))))
            print('==================')
            n_informative = nf // 10
            X, Y, coef_ = make_regression(n_samples=ns, n_features=nf,
                                          n_informative=n_informative,
                                          noise=0.1, n_targets=n_targets,
                                          coef=True)
            if Y.ndim == 1:
                Y = Y[:, None]
                coef_ = coef_[:, None]
            n_targets = Y.shape[1]

            X /= np.sqrt(np.sum(X ** 2, axis=0))  # Normalize data

            print("- benchmarking Lasso (%s)" % backend)
            X = np.asfortranarray(X)
            l1_reg = ns * alpha * l1_ratio
            l2_reg = ns * alpha * (1. - l1_ratio)
            if backend == "legacy":
                solver = cd_fast
            else:
                solver = cd_fast2
            W = np.zeros((n_targets, nf), dtype=X.dtype, order="F")
            tstart = time()
            if multi_task:
                n_iter = solver.enet_coordinate_descent_multi_task(
                    W, l1_reg, l2_reg, X, Y, max_iter, tol, rng)[-1]
                results.append(time() - tstart)
            elif precompute:
                assert n_targets == 1
                Gram = np.ascontiguousarray(np.dot(X.T, X))
                Cov = np.dot(X.T, Y)
                W = W[0]
                Cov = Cov[:, 0]
                Y = Y[:, 0]
                n_iter = solver.enet_coordinate_descent_gram(
                    W, l1_reg, l2_reg, Gram, Cov, Y, max_iter, tol, rng)[-1]
            else:
                assert n_targets == 1
                Y = Y[:, 0]
                W = W[0]
                n_iter = solver.enet_coordinate_descent(
                    W, l1_reg, l2_reg, X, Y, max_iter, tol, rng)[-1]
            results.append(time() - tstart)
    return results


if __name__ == '__main__':
    from sklearn.linear_model import Lasso, MultiTaskLasso
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_style("darkgrid")

    alpha = 0.1  # regularization parameter
    l1_ratio = 1.

    plt.figure('scikit-learn LASSO benchmark results')
    for i, precompute in enumerate([True, False]):
        plt.subplot("21%i" % (i + 1))
        if precompute:
            n_features = 100
            list_n_features = [n_features]
            list_n_samples = np.linspace(100, 100000, 5).astype(np.int)
            plt.title('precomputed Gram matrix, %d features, alpha=%s' % (
                n_features, alpha))
            plt.xlabel('number of samples')
        else:
            n_samples = 2000
            list_n_samples = [n_samples]
            list_n_features = np.linspace(500, 4000, 5).astype(np.int)
            plt.title('%d samples, alpha=%s' % (n_samples, alpha))
            plt.xlabel('number of features')
        plt.ylabel('Time (s)')
        for backend, color in zip(["legacy", "ninja"],
                                  ["r", "b"]):
            for multi_task, style in zip([True, False],
                                         ["-", "--"]):
                if multi_task and precompute:
                    continue
                if multi_task:
                    tag = "Multi"
                else:
                    tag = "Single"
                results = compute_bench(
                    alpha, list_n_samples, list_n_features, l1_ratio=l1_ratio,
                    multi_task=multi_task, backend=backend,
                    n_targets=10 if multi_task else 1, precompute=precompute)
                if precompute:
                    x = list_n_samples
                else:
                    x = list_n_features
                plt.plot(x, results, '%s%s' % (color, style), linewidth=2,
                         label='%s-task Lasso (%s)' % (tag, backend))
        plt.axis('tight')
        plt.legend(loc="best")
    plt.tight_layout()
    # plt.savefig("paper/figs/bench.png", dpi=200, bbox_inches="tight")
    plt.show()
