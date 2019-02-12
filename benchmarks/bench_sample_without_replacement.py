"""
Benchmarks for sampling without replacement of integer.

"""
import gc
import sys
import optparse
from datetime import datetime
import operator

import matplotlib.pyplot as plt
import numpy as np
import random

from sklearn.utils.random import sample_without_replacement


def compute_time(t_start, delta):
    mu_second = 0.0 + 10 ** 6  # number of microseconds in a second

    return delta.seconds + delta.microseconds / mu_second


def bench_sample(sampling, n_population, n_samples):
    gc.collect()
    # start time
    t_start = datetime.now()
    sampling(n_population, n_samples)
    delta = (datetime.now() - t_start)
    # stop time
    time = compute_time(t_start, delta)
    return time

if __name__ == "__main__":
    ###########################################################################
    # Option parser
    ###########################################################################
    op = optparse.OptionParser()
    op.add_option("--n-times",
                  dest="n_times", default=5, type=int,
                  help="Benchmark results are average over n_times experiments")

    op.add_option("--n-population",
                  dest="n_population", default=100000, type=int,
                  help="Size of the population to sample from.")

    op.add_option("--n-step",
                  dest="n_steps", default=5, type=int,
                  help="Number of step interval between 0 and n_population.")

    default_algorithms = "custom-tracking-selection,custom-auto," \
                         "custom-reservoir-sampling,custom-pool,"\
                         "python-core-sample,numpy-permutation"

    op.add_option("--algorithm",
                  dest="selected_algorithm",
                  default=default_algorithms,
                  type=str,
                  help="Comma-separated list of transformer to benchmark. "
                       "Default: %default. \nAvailable: %default")

    # op.add_option("--random-seed",
    #               dest="random_seed", default=13, type=int,
    #               help="Seed used by the random number generators.")

    (opts, args) = op.parse_args()
    if len(args) > 0:
        op.error("this script takes no arguments.")
        sys.exit(1)

    selected_algorithm = opts.selected_algorithm.split(',')
    for key in selected_algorithm:
        if key not in default_algorithms.split(','):
            raise ValueError("Unknown sampling algorithm \"%s\" not in (%s)."
                             % (key, default_algorithms))

    ###########################################################################
    # List sampling algorithm
    ###########################################################################
    # We assume that sampling algorithm has the following signature:
    #   sample(n_population, n_sample)
    #
    sampling_algorithm = {}

    ###########################################################################
    # Set Python core input
    sampling_algorithm["python-core-sample"] = \
        lambda n_population, n_sample: \
        random.sample(range(n_population), n_sample)

    ###########################################################################
    # Set custom automatic method selection
    sampling_algorithm["custom-auto"] = \
        lambda n_population, n_samples, random_state=None: \
        sample_without_replacement(n_population, n_samples, method="auto",
                                   random_state=random_state)

    ###########################################################################
    # Set custom tracking based method
    sampling_algorithm["custom-tracking-selection"] = \
        lambda n_population, n_samples, random_state=None: \
        sample_without_replacement(n_population,
                                   n_samples,
                                   method="tracking_selection",
                                   random_state=random_state)

    ###########################################################################
    # Set custom reservoir based method
    sampling_algorithm["custom-reservoir-sampling"] = \
        lambda n_population, n_samples, random_state=None: \
        sample_without_replacement(n_population,
                                   n_samples,
                                   method="reservoir_sampling",
                                   random_state=random_state)

    ###########################################################################
    # Set custom reservoir based method
    sampling_algorithm["custom-pool"] = \
        lambda n_population, n_samples, random_state=None: \
        sample_without_replacement(n_population,
                                   n_samples,
                                   method="pool",
                                   random_state=random_state)

    ###########################################################################
    # Numpy permutation based
    sampling_algorithm["numpy-permutation"] = \
        lambda n_population, n_sample: \
        np.random.permutation(n_population)[:n_sample]

    ###########################################################################
    # Remove unspecified algorithm
    sampling_algorithm = {key: value
                          for key, value in sampling_algorithm.items()
                          if key in selected_algorithm}

    ###########################################################################
    # Perform benchmark
    ###########################################################################
    time = {}
    n_samples = np.linspace(start=0, stop=opts.n_population,
        num=opts.n_steps).astype(np.int)

    ratio = n_samples / opts.n_population

    print('Benchmarks')
    print("===========================")

    for name in sorted(sampling_algorithm):
        print("Perform benchmarks for %s..." % name, end="")
        time[name] = np.zeros(shape=(opts.n_steps, opts.n_times))

        for step in range(opts.n_steps):
            for it in range(opts.n_times):
                time[name][step, it] = bench_sample(sampling_algorithm[name],
                                                    opts.n_population,
                                                    n_samples[step])

        print("done")

    print("Averaging results...", end="")
    for name in sampling_algorithm:
        time[name] = np.mean(time[name], axis=1)
    print("done\n")

    # Print results
    ###########################################################################
    print("Script arguments")
    print("===========================")
    arguments = vars(opts)
    print("%s \t | %s " % ("Arguments".ljust(16),
                           "Value".center(12),))
    print(25 * "-" + ("|" + "-" * 14) * 1)
    for key, value in arguments.items():
        print("%s \t | %s " % (str(key).ljust(16),
                               str(value).strip().center(12)))
    print("")

    print("Sampling algorithm performance:")
    print("===============================")
    print("Results are averaged over %s repetition(s)." % opts.n_times)
    print("")

    fig = plt.figure('scikit-learn sample w/o replacement benchmark results')
    plt.title("n_population = %s, n_times = %s" %
              (opts.n_population, opts.n_times))
    ax = fig.add_subplot(111)
    for name in sampling_algorithm:
        ax.plot(ratio, time[name], label=name)

    ax.set_xlabel('ratio of n_sample / n_population')
    ax.set_ylabel('Time (s)')
    ax.legend()

    # Sort legend labels
    handles, labels = ax.get_legend_handles_labels()
    hl = sorted(zip(handles, labels), key=operator.itemgetter(1))
    handles2, labels2 = zip(*hl)
    ax.legend(handles2, labels2, loc=0)

    plt.show()
