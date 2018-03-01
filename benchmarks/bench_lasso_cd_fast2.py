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
                  n_targets=10, max_iter=1000, tol=0, l1_ratio=1.,
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
                if precompute:
                    Gram = np.ascontiguousarray(np.dot(X.T, X))
                    Cov = np.dot(X.T, Y)
                    n_iter = solver.enet_coordinate_descent_multi_task_gram(
                        W, l1_reg, l2_reg, Gram, Cov, Y, max_iter, tol,
                        rng)[-1]
                else:
                    n_iter = solver.enet_coordinate_descent_multi_task(
                        W, l1_reg, l2_reg, X, Y, max_iter, tol, rng)[-1]
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
            if tol <= 0.:
                assert(n_iter == max_iter)
    return results


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_style("darkgrid")

    alpha = 0.01  # regularization parameter
    l1_ratio = 1.

    plt.figure('scikit-learn LASSO benchmark results')
    for i, precompute in enumerate([True, False]):
        if precompute:
            n_features = 50
            list_n_features = [n_features]
            list_n_samples = np.linspace(100, 1000000, 5).astype(np.int)
        else:
            n_samples = 200
            list_n_samples = [n_samples]
            list_n_features = np.linspace(500, 5000, 5).astype(np.int)
        for j, multi_task in enumerate([False, True]):
            plt.subplot2grid((2, 2), (i, j))
            if j == 0:
                plt.ylabel('Time (s)')
            if precompute:
                if multi_task:
                    plt.title("Multi-task 'Gram mode' (nf=%i)" % n_features)
                else:
                    plt.title("Single-task 'Gram mode' (nf=%i)" % n_features)
                plt.xlabel('number of samples')
            else:
                if multi_task:
                    plt.title("Multi-task 'Non-Gram mode' (ns=%i)" % n_samples)
                else:
                    plt.title("Single-task 'Non-Gram mode' (ns=%i)" % (
                        n_samples))
                plt.xlabel('number of features')
            for backend, color in zip(["legacy", "ninja"], ["r", "b"]):
                if multi_task and precompute and backend == "legacy":
                    continue
                results = compute_bench(
                    alpha, list_n_samples, list_n_features, l1_ratio=l1_ratio,
                    multi_task=multi_task, backend=backend,
                    n_targets=5 if multi_task else 1, precompute=precompute)
                if precompute:
                    x = list_n_samples
                else:
                    x = list_n_features
                plt.plot(x, results, "o-", c=color, linewidth=3,
                         label='Lasso (%s)' % backend)
    plt.legend(loc="best")
    plt.tight_layout()
    out_file = "bench.png"
    plt.savefig(out_file, dpi=200, bbox_inches="tight")
    print(out_file)
    plt.show()
