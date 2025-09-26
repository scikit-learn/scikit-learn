"""
=======================================================================
Benchmark: Comparing Python Wrappers for libsvm
=======================================================================

This script benchmarks the performance of three different Python wrappers for
the popular SVM implementation, libsvm:

1.  scikit-learn (sklearn.svm.SVC)
2.  The official SWIG-based bindings included with libsvm (svmutil)
3.  PyMVPA's wrapper (mvpa.clfs.svm)

Two benchmarks are performed:
----------------------------
1.  **Varying Samples:** The number of features is fixed, and we plot the
    execution time as the number of training samples increases.

2.  **Varying Dimensions:** The number of samples is fixed, and we plot the
    execution time as the number of features (dimensions) increases.

Requirements:
-------------
  * numpy
  * matplotlib
  * scikit-learn
  * libsvm (with its Python bindings)
  * pymvpa2
"""
import time
import numpy as np
import matplotlib.pyplot as plt

# --- Check for Dependencies ---
try:
    from sklearn.svm import SVC
except ImportError:
    print("Error: scikit-learn is not installed. Please run 'pip install scikit-learn'")
    exit()

try:
    import svmutil
except ImportError:
    print("Error: libsvm python bindings are not found.")
    print("Please install libsvm and ensure the 'svmutil.py' is in your PYTHONPATH.")
    exit()

try:
    from mvpa2.datasets import Dataset
    from mvpa2.clfs import svm as mvpa_svm
except ImportError:
    print("Error: PyMVPA is not installed. Please run 'pip install pymvpa2'")
    exit()


# ==============================================================================
# Benchmarking Functions
# ==============================================================================

def benchmark_sklearn(X: np.ndarray, y: np.ndarray) -> float:
    """Benchmarks scikit-learn's SVC."""
    t_start = time.perf_counter()
    clf = SVC(kernel='rbf', gamma='auto')
    clf.fit(X, y).predict(X)
    return time.perf_counter() - t_start

def benchmark_libsvm(X: np.ndarray, y: np.ndarray) -> float:
    """Benchmarks the official libsvm python bindings."""
    # libsvm requires lists instead of numpy arrays
    X_list = X.tolist()
    y_list = y.tolist()

    t_start = time.perf_counter()
    problem = svmutil.svm_problem(y_list, X_list)
    # Parameters for an RBF kernel C-SVC
    param = svmutil.svm_parameter('-s 0 -t 2 -q')
    model = svmutil.svm_train(problem, param)
    svmutil.svm_predict(y_list, X_list, model)
    return time.perf_counter() - t_start

def benchmark_pymvpa(X: np.ndarray, y: np.ndarray) -> float:
    """Benchmarks PyMVPA's libsvm wrapper."""
    t_start = time.perf_counter()
    dataset = Dataset(samples=X, labels=y)
    clf = mvpa_svm.RbfCSVMC(C=1.0)
    clf.train(dataset)
    clf.predict(X)
    return time.perf_counter() - t_start


# ==============================================================================
# Plotting Function
# ==============================================================================

def plot_results(ax, x_values, results: dict, title: str, xlabel: str):
    """
    Helper function to plot the benchmark results on a given axis.
    """
    ax.set_title(title, fontsize=14)
    ax.plot(x_values, results['PyMVPA'], 'g-o', label='PyMVPA')
    ax.plot(x_values, results['libsvm'], 'r-s', label='libsvm (Official)')
    ax.plot(x_values, results['scikit-learn'], 'b-^', label='scikit-learn')
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Time (seconds)')
    ax.legend()
    ax.grid(True)


# ==============================================================================
# Main Benchmark Drivers
# ==============================================================================

def run_benchmark_by_samples(ax, n_iterations=5, start_samples=200, step=100, dim=200):
    """Runs benchmark with a varying number of samples."""
    print("\n--- Starting Benchmark 1: Varying Number of Samples ---")
    results = {'scikit-learn': [], 'libsvm': [], 'PyMVPA': []}
    sample_sizes = []

    for i in range(n_iterations):
        n_samples = start_samples + i * step
        sample_sizes.append(n_samples)
        print(f"Iteration {i+1}/{n_iterations}: Processing {n_samples} samples...")

        X = np.random.randn(n_samples, dim)
        y = (np.random.randn(n_samples) * 10).astype(int)

        results['scikit-learn'].append(benchmark_sklearn(X, y))
        results['libsvm'].append(benchmark_libsvm(X, y))
        results['PyMVPA'].append(benchmark_pymvpa(X, y))

    plot_results(ax, sample_sizes, results,
                 'SVM Performance vs. Number of Samples',
                 'Number of Samples to Train and Classify')

def run_benchmark_by_dimensions(ax, n_iterations=10, n_samples=100, start_dim=100, step=500):
    """Runs benchmark with a varying number of dimensions."""
    print("\n--- Starting Benchmark 2: Varying Number of Dimensions ---")
    print("Warning: This may take a long time.")
    results = {'scikit-learn': [], 'libsvm': [], 'PyMVPA': []}
    dimensions = []

    for i in range(n_iterations):
        dim = start_dim + i * step
        dimensions.append(dim)
        print(f"Iteration {i+1}/{n_iterations}: Processing {dim} dimensions...")

        X = np.random.randn(n_samples, dim)
        y = (np.random.randn(n_samples) * 10).astype(int)

        results['scikit-learn'].append(benchmark_sklearn(X, y))
        results['libsvm'].append(benchmark_libsvm(X, y))
        results['PyMVPA'].append(benchmark_pymvpa(X, y))

    plot_results(ax, dimensions, results,
                 'SVM Performance in High Dimensional Spaces',
                 'Number of Dimensions')


if __name__ == '__main__':
    print(__doc__)

    # Create a figure with two subplots (one for each benchmark)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))

    # Run the benchmarks
    run_benchmark_by_samples(ax1)
    run_benchmark_by_dimensions(ax2)

    # Show the final plot
    fig.tight_layout(pad=3.0)
    plt.show()