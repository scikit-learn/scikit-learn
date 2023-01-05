#!/usr/bin/env python
"""
A comparison of multilabel target formats and metrics over them
"""

from timeit import timeit
from functools import partial
import itertools
import argparse
import sys

import matplotlib.pyplot as plt
import scipy.sparse as sp
import numpy as np

from sklearn.datasets import make_multilabel_classification
from sklearn.metrics import (
    f1_score,
    accuracy_score,
    hamming_loss,
    jaccard_similarity_score,
)
from sklearn.utils._testing import ignore_warnings


METRICS = {
    "f1": partial(f1_score, average="micro"),
    "f1-by-sample": partial(f1_score, average="samples"),
    "accuracy": accuracy_score,
    "hamming": hamming_loss,
    "jaccard": jaccard_similarity_score,
}

FORMATS = {
    "sequences": lambda y: [list(np.flatnonzero(s)) for s in y],
    "dense": lambda y: y,
    "csr": sp.csr_matrix,
    "csc": sp.csc_matrix,
}


@ignore_warnings
def benchmark(
    metrics=tuple(v for k, v in sorted(METRICS.items())),
    formats=tuple(v for k, v in sorted(FORMATS.items())),
    samples=1000,
    classes=4,
    density=0.2,
    n_times=5,
):
    """Times metric calculations for a number of inputs

    Parameters
    ----------
    metrics : array-like of callables (1d or 0d)
        The metric functions to time.

    formats : array-like of callables (1d or 0d)
        These may transform a dense indicator matrix into multilabel
        representation.

    samples : array-like of ints (1d or 0d)
        The number of samples to generate as input.

    classes : array-like of ints (1d or 0d)
        The number of classes in the input.

    density : array-like of ints (1d or 0d)
        The density of positive labels in the input.

    n_times : int
        Time calling the metric n_times times.

    Returns
    -------
    array of floats shaped like (metrics, formats, samples, classes, density)
        Time in seconds.
    """
    metrics = np.atleast_1d(metrics)
    samples = np.atleast_1d(samples)
    classes = np.atleast_1d(classes)
    density = np.atleast_1d(density)
    formats = np.atleast_1d(formats)
    out = np.zeros(
        (len(metrics), len(formats), len(samples), len(classes), len(density)),
        dtype=float,
    )
    it = itertools.product(samples, classes, density)
    for i, (s, c, d) in enumerate(it):
        _, y_true = make_multilabel_classification(
            n_samples=s, n_features=1, n_classes=c, n_labels=d * c, random_state=42
        )
        _, y_pred = make_multilabel_classification(
            n_samples=s, n_features=1, n_classes=c, n_labels=d * c, random_state=84
        )
        for j, f in enumerate(formats):
            f_true = f(y_true)
            f_pred = f(y_pred)
            for k, metric in enumerate(metrics):
                t = timeit(partial(metric, f_true, f_pred), number=n_times)

                out[k, j].flat[i] = t
    return out


def _tabulate(results, metrics, formats):
    """Prints results by metric and format

    Uses the last ([-1]) value of other fields
    """
    column_width = max(max(len(k) for k in formats) + 1, 8)
    first_width = max(len(k) for k in metrics)
    head_fmt = "{:<{fw}s}" + "{:>{cw}s}" * len(formats)
    row_fmt = "{:<{fw}s}" + "{:>{cw}.3f}" * len(formats)
    print(head_fmt.format("Metric", *formats, cw=column_width, fw=first_width))
    for metric, row in zip(metrics, results[:, :, -1, -1, -1]):
        print(row_fmt.format(metric, *row, cw=column_width, fw=first_width))


def _plot(
    results,
    metrics,
    formats,
    title,
    x_ticks,
    x_label,
    format_markers=("x", "|", "o", "+"),
    metric_colors=("c", "m", "y", "k", "g", "r", "b"),
):
    """
    Plot the results by metric, format and some other variable given by
    x_label
    """
    fig = plt.figure("scikit-learn multilabel metrics benchmarks")
    plt.title(title)
    ax = fig.add_subplot(111)
    for i, metric in enumerate(metrics):
        for j, format in enumerate(formats):
            ax.plot(
                x_ticks,
                results[i, j].flat,
                label="{}, {}".format(metric, format),
                marker=format_markers[j],
                color=metric_colors[i % len(metric_colors)],
            )
    ax.set_xlabel(x_label)
    ax.set_ylabel("Time (s)")
    ax.legend()
    plt.show()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "metrics",
        nargs="*",
        default=sorted(METRICS),
        help="Specifies metrics to benchmark, defaults to all. Choices are: {}".format(
            sorted(METRICS)
        ),
    )
    ap.add_argument(
        "--formats",
        nargs="+",
        choices=sorted(FORMATS),
        help="Specifies multilabel formats to benchmark (defaults to all).",
    )
    ap.add_argument(
        "--samples", type=int, default=1000, help="The number of samples to generate"
    )
    ap.add_argument("--classes", type=int, default=10, help="The number of classes")
    ap.add_argument(
        "--density",
        type=float,
        default=0.2,
        help="The average density of labels per sample",
    )
    ap.add_argument(
        "--plot",
        choices=["classes", "density", "samples"],
        default=None,
        help=(
            "Plot time with respect to this parameter varying up to the specified value"
        ),
    )
    ap.add_argument(
        "--n-steps", default=10, type=int, help="Plot this many points for each metric"
    )
    ap.add_argument(
        "--n-times", default=5, type=int, help="Time performance over n_times trials"
    )
    args = ap.parse_args()

    if args.plot is not None:
        max_val = getattr(args, args.plot)
        if args.plot in ("classes", "samples"):
            min_val = 2
        else:
            min_val = 0
        steps = np.linspace(min_val, max_val, num=args.n_steps + 1)[1:]
        if args.plot in ("classes", "samples"):
            steps = np.unique(np.round(steps).astype(int))
        setattr(args, args.plot, steps)

    if args.metrics is None:
        args.metrics = sorted(METRICS)
    if args.formats is None:
        args.formats = sorted(FORMATS)

    results = benchmark(
        [METRICS[k] for k in args.metrics],
        [FORMATS[k] for k in args.formats],
        args.samples,
        args.classes,
        args.density,
        args.n_times,
    )

    _tabulate(results, args.metrics, args.formats)

    if args.plot is not None:
        print("Displaying plot", file=sys.stderr)
        title = "Multilabel metrics with %s" % ", ".join(
            "{0}={1}".format(field, getattr(args, field))
            for field in ["samples", "classes", "density"]
            if args.plot != field
        )
        _plot(results, args.metrics, args.formats, title, steps, args.plot)
