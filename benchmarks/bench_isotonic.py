"""
Benchmarks of isotonic regression performance.

We generate a synthetic dataset of size 10^n, for n in [min, max], and
examine the time taken to run isotonic regression over the dataset.

The timings are then output to stdout, or visualized on a log-log scale
with matplotlib.

This allows the scaling of the algorithm with the problem size to be
visualized and understood.
"""
import argparse
import gc
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import expit

from sklearn.isotonic import isotonic_regression


def generate_perturbed_logarithm_dataset(size):
    return np.random.randint(-50, 50, size=size) + 50.0 * np.log(1 + np.arange(size))


def generate_logistic_dataset(size):
    X = np.sort(np.random.normal(size=size))
    return np.random.random(size=size) < expit(X)


def generate_pathological_dataset(size):
    # Triggers O(n^2) complexity on the original implementation.
    return np.r_[
        np.arange(size), np.arange(-(size - 1), size), np.arange(-(size - 1), 1)
    ]


DATASET_GENERATORS = {
    "perturbed_logarithm": generate_perturbed_logarithm_dataset,
    "logistic": generate_logistic_dataset,
    "pathological": generate_pathological_dataset,
}


def bench_isotonic_regression(Y):
    """
    Runs a single iteration of isotonic regression on the input data,
    and reports the total time taken (in seconds).
    """
    gc.collect()

    tstart = datetime.now()
    isotonic_regression(Y)
    return (datetime.now() - tstart).total_seconds()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Isotonic Regression benchmark tool")
    parser.add_argument("--seed", type=int, help="RNG seed")
    parser.add_argument(
        "--iterations",
        type=int,
        required=True,
        help="Number of iterations to average timings over for each problem size",
    )
    parser.add_argument(
        "--log_min_problem_size",
        type=int,
        required=True,
        help="Base 10 logarithm of the minimum problem size",
    )
    parser.add_argument(
        "--log_max_problem_size",
        type=int,
        required=True,
        help="Base 10 logarithm of the maximum problem size",
    )
    parser.add_argument(
        "--show_plot", action="store_true", help="Plot timing output with matplotlib"
    )
    parser.add_argument("--dataset", choices=DATASET_GENERATORS.keys(), required=True)

    args = parser.parse_args()

    np.random.seed(args.seed)

    timings = []
    for exponent in range(args.log_min_problem_size, args.log_max_problem_size):
        n = 10**exponent
        Y = DATASET_GENERATORS[args.dataset](n)
        time_per_iteration = [
            bench_isotonic_regression(Y) for i in range(args.iterations)
        ]
        timing = (n, np.mean(time_per_iteration))
        timings.append(timing)

        # If we're not plotting, dump the timing to stdout
        if not args.show_plot:
            print(n, np.mean(time_per_iteration))

    if args.show_plot:
        plt.plot(*zip(*timings))
        plt.title("Average time taken running isotonic regression")
        plt.xlabel("Number of observations")
        plt.ylabel("Time (s)")
        plt.axis("tight")
        plt.loglog()
        plt.show()
