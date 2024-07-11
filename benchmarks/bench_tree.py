"""
To run this, you'll need to have installed.

  * scikit-learn

Does two benchmarks

First, we fix a training set, increase the number of
samples to classify and plot number of classified samples as a
function of time.

In the second benchmark, we increase the number of dimensions of the
training set, classify a sample and plot the time taken as a function
of the number of dimensions.
"""

import gc
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np

# to store the results
scikit_classifier_results = []
scikit_regressor_results = []
meta_list = []

mu_second = 0.0 + 10**6  # number of microseconds in a second


def bench_scikit_tree_classifier(X, Y):
    """Benchmark with scikit-learn decision tree classifier"""

    from sklearn.tree import DecisionTreeClassifier

    gc.collect()

    # start time
    tstart = datetime.now()
    clf = DecisionTreeClassifier()
    clf.fit(X, Y).predict(X)
    delta = datetime.now() - tstart
    # stop time

    scikit_classifier_results.append(delta.seconds + delta.microseconds / mu_second)


def bench_scikit_tree_regressor(X, Y):
    """Benchmark with scikit-learn decision tree regressor"""

    from sklearn.tree import DecisionTreeRegressor

    gc.collect()

    # start time
    tstart = datetime.now()
    clf = DecisionTreeRegressor()
    clf.fit(X, Y).predict(X)
    delta = datetime.now() - tstart
    # stop time

    scikit_regressor_results.append(delta.seconds + delta.microseconds / mu_second)


def plot():
    import pandas as pd
    import seaborn as sns
    
    # Function to load and process results
    def load_results(branch, result_type):
        results = np.load(f"scikit_tree_{result_type}_results_{branch}.npz")
        classifier_results = results['scikit_classifier_results']
        regressor_results = results['scikit_regressor_results']
        meta_list = results['meta_list']
        return classifier_results, regressor_results, meta_list

    # Load the results for both branches
    main_clf_samples, main_reg_samples, main_meta_samples = load_results('main', 'samples')
    pr_clf_samples, pr_reg_samples, pr_meta_samples = load_results('pr_no_fusedtype', 'samples')
    main_clf_dims, main_reg_dims, main_meta_dims = load_results('main', 'dims')
    pr_clf_dims, pr_reg_dims, pr_meta_dims = load_results('pr_no_fusedtype', 'dims')

    # Ensure meta_lists are arrays
    main_meta_samples = np.array(main_meta_samples)
    pr_meta_samples = np.array(pr_meta_samples)
    main_meta_dims = np.array(main_meta_dims)
    pr_meta_dims = np.array(pr_meta_dims)

    # Check shapes of loaded arrays
    print(f"main_clf_samples shape: {main_clf_samples.shape}")
    print(f"main_reg_samples shape: {main_reg_samples.shape}")
    print(f"pr_clf_samples shape: {pr_clf_samples.shape}")
    print(f"pr_reg_samples shape: {pr_reg_samples.shape}")
    print(f"main_meta_samples shape: {main_meta_samples.shape}")
    print(f"pr_meta_samples shape: {pr_meta_samples.shape}")

    print(f"main_clf_dims shape: {main_clf_dims.shape}")
    print(f"main_reg_dims shape: {main_reg_dims.shape}")
    print(f"pr_clf_dims shape: {pr_clf_dims.shape}")
    print(f"pr_reg_dims shape: {pr_reg_dims.shape}")
    print(f"main_meta_dims shape: {main_meta_dims.shape}")
    print(f"pr_meta_dims shape: {pr_meta_dims.shape}")

    # Create DataFrames for samples and dimensions
    samples_df = pd.DataFrame({
        'n_samples': main_meta_samples,
        'classifier_main': main_clf_samples.flatten(),
        'regressor_main': main_reg_samples.flatten(),
        'classifier_pr': pr_clf_samples.flatten(),
        'regressor_pr': pr_reg_samples.flatten()
    })

    dims_df = pd.DataFrame({
        'n_dims': main_meta_dims,
        'classifier_main': main_clf_dims.flatten(),
        'regressor_main': main_reg_dims.flatten(),
        'classifier_pr': pr_clf_dims.flatten(),
        'regressor_pr': pr_reg_dims.flatten()
    })

    # Melt the DataFrames for easier plotting
    samples_df = samples_df.melt(id_vars='n_samples', var_name='method', value_name='time')
    dims_df = dims_df.melt(id_vars='n_dims', var_name='method', value_name='time')

    # Plotting with seaborn
    plt.figure(figsize=(14, 10))

    # Plot for varying number of samples
    plt.subplot(211)
    sns.lineplot(data=samples_df, x='n_samples', y='time', hue='method', ci='sd')
    plt.title("Learning with varying number of samples")
    plt.xlabel("Number of Samples")
    plt.ylabel("Time (s)")

    # Plot for varying number of dimensions
    plt.subplot(212)
    sns.lineplot(data=dims_df, x='n_dims', y='time', hue='method', ci='sd')
    plt.title("Learning in high dimensional spaces")
    plt.xlabel("Number of Dimensions")
    plt.ylabel("Time (s)")

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    plot()
    assert False
    print("============================================")
    print("Warning: this is going to take a looong time")
    print("============================================")
    branch = 'pr_no_fusedtype'
    n_repeat = 5
    n = 10
    step = 10000
    n_samples = 10000
    dim = 10
    n_classes = 10
    for _ in range(n_repeat):
        for i in range(n):
            print("============================================")
            print("Entering iteration %s of %s" % (i, n))
            print("============================================")
            n_samples += step
            X = np.random.randn(n_samples, dim)
            Y = np.random.randint(0, n_classes, (n_samples,))
            bench_scikit_tree_classifier(X, Y)
            Y = np.random.randn(n_samples)
            bench_scikit_tree_regressor(X, Y)
            meta_list.append(n_samples)

    # xx = range(0, n * step, step)
    # plt.figure("scikit-learn tree benchmark results")
    # plt.subplot(211)
    # plt.title("Learning with varying number of samples")
    # plt.plot(xx, scikit_classifier_results, "g-", label="classification")
    # plt.plot(xx, scikit_regressor_results, "r-", label="regression")
    # plt.legend(loc="upper left")
    # plt.xlabel("number of samples")
    # plt.ylabel("Time (s)")

    # XXX: remove
    # save the results to a file
    np.savez(
        f"scikit_tree_samples_results_{branch}.npz",
        scikit_classifier_results=scikit_classifier_results,
        scikit_regressor_results=scikit_regressor_results,
        meta_list=meta_list,
    )

    scikit_classifier_results = []
    scikit_regressor_results = []
    meta_list = []
    n = 10
    step = 500
    start_dim = 500
    n_classes = 10

    dim = start_dim
    for _ in range(n_repeat):
        for i in range(0, n):
            print("============================================")
            print("Entering iteration %s of %s" % (i, n))
            print("============================================")
            dim += step
            X = np.random.randn(100, dim)
            Y = np.random.randint(0, n_classes, (100,))
            bench_scikit_tree_classifier(X, Y)
            Y = np.random.randn(100)
            bench_scikit_tree_regressor(X, Y)
            meta_list.append(dim)

    # xx = np.arange(start_dim, start_dim + n * step, step)
    # plt.subplot(212)
    # plt.title("Learning in high dimensional spaces")
    # plt.plot(xx, scikit_classifier_results, "g-", label="classification")
    # plt.plot(xx, scikit_regressor_results, "r-", label="regression")
    # plt.legend(loc="upper left")
    # plt.xlabel("number of dimensions")
    # plt.ylabel("Time (s)")
    # plt.axis("tight")
    # plt.show()

    # XXX: remove
    # save the results to a file
    np.savez(
        f"scikit_tree_dims_results_{branch}.npz",
        scikit_classifier_results=scikit_classifier_results,
        scikit_regressor_results=scikit_regressor_results,
        meta_list=meta_list,
    )
