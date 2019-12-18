"""
=================================================
Label Propagation MNIST: Comparing Sparse Kernels
=================================================

This example compares the runtime and performance of two sparse kernels for
semisupervised learning on the MNIST digit dataset.

The MNIST dataset consists of 28x28 pixel grayscale images. Here, we will use a
subset of 10K images, reserving a fraction of these for testing.  We will
compare the performance and runtime of two sparse kernels, across a range of
low-supervision scenarios.

In each scenario, we will run each model multiple times, to increase our
confidence in the comparison between kernels.

The models will be evaluated for their accuracy at spreading labels during
training ("transductive learning"), as well as spreading labels to unseen
points at test time ("inductive learning").

The first kernel option produces a binary k-Nearest Neighbors adjacency matrix.
The second produces a kernel which is also k-sparse, but contains the same
weights as used in an RBF kernel.

Notice that the performance of the sparse-RBF kernel is very sensitive to
parameters; the parameters used here were found by a quick manual search, so
the model can likely be improved with further optimization, and using this
kernel effectively on a new dataset requires hyperparameter tuning.
"""
import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint
from sklearn.datasets import fetch_openml
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.semi_supervised import LabelSpreading
from torchvision.datasets import CIFAR10
import time


def run_grid_search(X_train, X_test, y_train, y_test):
    """
    We perform a grid search to optimize parameters for sparse-rbf
    kernel.  For this purpose, we use a smaller subset of the data.  In all
    cases, we simply use max_iter=100 Notice that we are searching over the
    inductive accuracy (accuracy on the test set) rather than the
    transductive accuracy (on the masked training examples).  This keeps
    things a bit simpler, though we could customize the score function and
    the `WrapLabelSpreading` class further to also hyperparameter search over
    the transductive accuracy.
    """
    # In order to use GridSearchCV, we will use a thin wrapper class
    # that masks a subset of our training labels.
    sparse_rbf_model = GridSearchCV(
            WrapLabelSpreading(kernel='sparse-rbf', supervision_fraction=0.05),
            param_grid={
                'n_jobs': [-1],
                'max_iter': [100],
                'alpha': np.linspace(0.01, 0.99, 10),
                'gamma': np.logspace(-8, -4, 10),
                'n_neighbors': list(range(6, 30, 2))},
            cv=3)

    sparse_rbf_model.fit(np.vstack((X_train, X_test)),
                         np.concatenate((y_train, y_test)))
    sparse_rbf_params = sparse_rbf_model.best_params_
    print(f"Optimal parameters for sparse RBF kernel: {sparse_rbf_params}")

    knn_model = GridSearchCV(
            WrapLabelSpreading(kernel='knn', supervision_fraction=0.05),
            param_grid={
                'n_jobs': [-1],
                'max_iter': [100],
                'alpha': np.linspace(0.01, 0.99, 10),
                'n_neighbors': list(range(6, 30, 2))},
            cv=3)

    knn_model.fit(np.vstack((X_train, X_test)),
                  np.concatenate((y_train, y_test)))
    knn_params = knn_model.best_params_
    print(f"Optimal parameters for knn kernel: {knn_params}")
    return sparse_rbf_params, knn_params


def run_comparison(X_train, X_test, y_train, y_test,
                   sparse_rbf_params, knn_params):
    print("Begin comparison...")
    supervision_fractions = [0.001, 0.003, 0.005, 0.01, 0.03, 0.05, 0.1]
    results = {
            'transduction': {
                'knn': {'avg': [], 'std': []},
                'sparse-rbf': {'avg': [], 'std': []}},
            'induction': {
                'knn': {'avg': [], 'std': []},
                'sparse-rbf': {'avg': [], 'std': []}},
            'runtimes': {
                'knn': {'avg': [], 'std': []},
                'sparse-rbf': {'avg': [], 'std': []}}
    }

    rng = np.random.RandomState(0)
    for supervision_fraction in supervision_fractions:
        for kernel_name, params in zip(['knn', 'sparse-rbf'],
                                       [knn_params, sparse_rbf_params]):
            model = LabelSpreading(kernel=kernel_name, **params)
            print("="*80)
            print(f"Kernel: {kernel_name}, " +
                  f"Supervision fraction: {supervision_fraction}")

            # Repeat each scenario several times to collect rough statistics
            t_accs = []
            i_accs = []
            runtimes = []
            for _ in range(5):
                n_train = len(y_train)
                n_labeled = int(supervision_fraction * n_train)
                indices = np.arange(n_train)
                rng.shuffle(indices)
                unlabeled_set = indices[n_labeled:]

                y_masked = np.copy(y_train)
                y_masked[unlabeled_set] = -1

                t0 = time.time()
                model.fit(X_train, y_masked)

                predicted_labels = model.transduction_[unlabeled_set]
                true_labels = y_train[unlabeled_set]
                t_acc = (np.sum(predicted_labels == true_labels) /
                         len(unlabeled_set))
                t_accs.append(t_acc)
                i_accs.append(model.score(X_test, y_test))
                t1 = time.time()
                runtimes.append(t1-t0)

            results['transduction'][kernel_name]['avg'].append(np.mean(t_accs))
            results['transduction'][kernel_name]['std'].append(np.std(t_accs))
            results['induction'][kernel_name]['avg'].append(np.mean(i_accs))
            results['induction'][kernel_name]['std'].append(np.std(i_accs))
            results['runtimes'][kernel_name]['avg'].append(np.mean(runtimes))
            results['runtimes'][kernel_name]['std'].append(np.std(runtimes))

    print("="*80)
    print(f"supervision_fractions: {supervision_fractions}")
    pprint(results)
    return supervision_fractions, results


def plot_results(supervision_fractions, results):
    fig, ax = plt.subplots(3, 1, figsize=(16, 9))
    for i, (label, ylabel) in enumerate(zip(
                ['induction', 'transduction', 'runtimes'],
                ['% Accuracy', '% Accuracy', 'Duration (s)'])):

        S_avg = results[label]['sparse-rbf']['avg']
        S_std = results[label]['sparse-rbf']['std']

        K_avg = results[label]['knn']['avg']
        K_std = results[label]['knn']['std']

        ax[i].scatter(supervision_fractions, S_avg, c='b', label='sparse-rbf')
        ax[i].scatter(supervision_fractions, K_avg, c='r', label='knn')
        ax[i].set_xscale('log')
        ax[i].set_xlim([8e-4, 1.3e-1])
        ax[i].set_title(f'{label.capitalize()}')
        ax[i].set_xlabel('Supervision Fraction')
        ax[i].fill_between(supervision_fractions,
                           [a - b for a, b in zip(S_avg, S_std)],
                           [a + b for a, b in zip(S_avg, S_std)],
                           facecolor='b', alpha=0.2)
        ax[i].fill_between(supervision_fractions,
                           [a + b for a, b in zip(K_avg, K_std)],
                           [a - b for a, b in zip(K_avg, K_std)],
                           facecolor='r', alpha=0.2)

    plt.tight_layout()
    plt.savefig('sparse_kernel_comparison.png')


class WrapLabelSpreading(LabelSpreading):
    """
    In order to perform a grid search over this semi-supervised model,
    we need to provide a wrapper that masks a subset of the data before
    `fit` is called.
    """
    def __init__(self, supervision_fraction, kernel='sparse-rbf', gamma=20,
                 n_neighbors=7, alpha=0.2, max_iter=30, tol=1e-3, n_jobs=None):

        self.supervision_fraction = supervision_fraction

        super().__init__(kernel=kernel, gamma=gamma,
                         n_neighbors=n_neighbors, alpha=alpha,
                         max_iter=max_iter, tol=tol, n_jobs=n_jobs)

    def fit(self, X, y):
        # mask a random subset of labels, based on self.supervision_fraction
        n_total = len(y)
        n_labeled = int(self.supervision_fraction * n_total)

        indices = np.arange(n_total)
        np.random.seed(0)
        np.random.shuffle(indices)
        unlabeled_subset = indices[n_labeled:]

        y[unlabeled_subset] = -1

        super().fit(X, y)
        return self


if __name__ == '__main__':
    # Choose the dataset
    dataset = 'mnist'
    # Set the fraction of data to use for hyperparam tuning
    hyperp_tune_fraction = 0.1
    # Set the fraction of data to use for the final comparison
    compare_fraction = 0.1
    # Set this flag to run the grid search, which takes several hours
    do_grid_search = False

    if dataset == 'mnist':
        X, y = fetch_openml('mnist_784', version=1, return_X_y=True)
        y = y.astype(int)

        X_train, X_test, y_train, y_test = train_test_split(
                X, y, test_size=10000, random_state=0)

    elif dataset == 'cifar10':
        cifar10 = CIFAR10('.', download=True, train=True)
        X_train = cifar10.data.reshape(-1, 3*32*32)
        y_train = np.array(cifar10.targets)

        cifar10 = CIFAR10('.', download=True, train=False)
        X_test = cifar10.data.reshape(-1, 3*32*32)
        y_test = np.array(cifar10.targets)
    else:
        raise ValueError(f"dataset {dataset} not supported")

    print("Full dataset sizes: " +
          f"\nX_train {X_train.shape}" +
          f"\nX_test {X_test.shape}" +
          f"\ny_train {y_train.shape}" +
          f"\ny_test {y_test.shape}")

    if do_grid_search:
        # First, we use a small subset of the data to tune hyperparameters
        tr_tune = int(hyperp_tune_fraction * len(y_train))
        te_tune = int(hyperp_tune_fraction * len(y_test))
        print("# items for hyperparam tuning:" +
              f"train: {tr_tune}, test: {te_tune}")
        sparse_rbf_params, knn_params = run_grid_search(
                X_train[:tr_tune, :],
                X_test[:te_tune, :],
                y_train[:tr_tune],
                y_test[:te_tune])
    else:
        # Values found from running grid search previously on MNIST
        tr_tune = 0
        te_tune = 0
        sparse_rbf_params = {
                'alpha': 0.663,
                'gamma': 2.154e-7,
                'max_iter': 100,
                'n_jobs': -1,
                'n_neighbors': 20}
        knn_params = {
                'alpha': 0.772,
                'max_iter': 100,
                'n_jobs': -1,
                'n_neighbors': 6}

    # Skip the items used for hyperparam tuning
    tr_comp = int(compare_fraction * len(y_train)) + tr_tune
    te_comp = int(compare_fraction * len(y_test)) + te_tune
    supervision_fractions, results = run_comparison(
            X_train[tr_tune:tr_comp, :], X_test[te_tune:te_comp, :],
            y_train[tr_tune:tr_comp],    y_test[te_tune:te_comp],
            sparse_rbf_params=sparse_rbf_params,
            knn_params=knn_params)

    plot_results(supervision_fractions, results)
