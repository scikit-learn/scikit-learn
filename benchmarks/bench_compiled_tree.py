from __future__ import print_function

from collections import namedtuple
from datetime import datetime
from functools import partial
from sklearn import tree, ensemble, datasets
from sklearn.tree.tree import DTYPE
from sklearn.utils import array2d
from sklearn.utils.bench import total_seconds
import argparse
import gc
import matplotlib.pyplot as plt
import numpy as np
import json


def with_depth(cls, max_depth):
    return partial(cls, max_depth=max_depth)


def with_best_first(cls, max_leaf_nodes):
    return partial(cls, max_leaf_nodes=max_leaf_nodes)


def uniform_dataset(args):
    X = np.random.random(size=(args.num_examples, args.num_features))
    y = np.random.choice([-1, 1], size=args.num_examples)
    return (X, y)

DATASETS = {
    "uniform": uniform_dataset,
    "hastie": lambda args: datasets.make_hastie_10_2(
        n_samples=args.num_examples),
    "friedman1": lambda args: datasets.make_friedman1(
        n_samples=args.num_examples, n_features=args.num_features),
    "friedman2": lambda args: datasets.make_friedman2(
        n_samples=args.num_examples, noise=args.noise),
    "friedman3": lambda args: datasets.make_friedman3(
        n_samples=args.num_examples, noise=args.noise),
    "make_regression": lambda args: datasets.make_regression(
        n_samples=args.num_examples,
        n_features=args.num_features,
        n_informative=args.num_informative)
}

ENSEMBLE_REGRESSORS = [
    ("GB-D1", with_depth(ensemble.GradientBoostingRegressor, 1)),
    ("GB-D3",  with_depth(ensemble.GradientBoostingRegressor, 3)),
    ("GB-B10", with_best_first(ensemble.GradientBoostingRegressor, 10)),
    ("RF-D1", with_depth(ensemble.RandomForestRegressor, 1)),
    ("RF-D3", with_depth(ensemble.RandomForestRegressor, 3)),
    ("RF-D5", with_depth(ensemble.RandomForestRegressor, 5)),
]

RelativeTiming = namedtuple('RelativeTiming', ['compiled', 'normal'])


def benchmark(args, predictor, X):
    total = 0.0
    for _ in range(args.iterations):
        gc.collect()
        start = datetime.now()
        predictor.predict(X)
        end = datetime.now()
        total += total_seconds(end - start)
    return total / args.iterations


def run_ensemble(args, name, cls, X, y):
    print(name)
    timing = lambda p: benchmark(args, p, X)

    def run(n_estimators):
        clf = cls(n_estimators=n_estimators)
        clf.fit(X, y)
        compiled = tree.CompiledRegressionPredictor(clf)
        relative_timing = RelativeTiming(
            compiled=timing(compiled), normal=timing(clf))
        print(n_estimators, relative_timing)
        return relative_timing

    return [(int(n), run(int(n)))
            for n in np.linspace(args.min_estimators,
                                 args.max_estimators,
                                 args.num_estimator_values)]


def titles(args):
    return (
        "Relative Speed of Compiled Evaluation", json.dumps(vars(args)))


def plot(args, timings):
    for name, cls_timings in timings:
        xs, relative_timings = zip(*cls_timings)
        ys = [r.normal / r.compiled for r in relative_timings]
        plt.plot(xs, ys, '-o', label=name)
        plt.hlines(1.0, np.min(xs), np.max(xs), 'k')

    plt.xlabel('Number of weak learners')
    plt.ylabel('Relative speedup')
    plt.axis('tight')
    plt.legend()
    plt.gca().set_ylim(bottom=0)
    title, suptitle = titles(args)
    plt.title(title)
    plt.suptitle(suptitle, fontsize=3)
    filename = "timings{0}.png".format(hash(str(args)))
    plt.savefig(filename, dpi=300)


def run_simulation(args):
    X, y = DATASETS[args.dataset](args)
    X = array2d(X, dtype=DTYPE)
    timings = [(name, run_ensemble(args, name, cls, X, y))
               for name, cls in ENSEMBLE_REGRESSORS]
    plot(args, timings)


def build_parser():
    parser = argparse.ArgumentParser(description="Timing ensemble evaluation")
    parser.add_argument('--iterations', type=int, required=True)
    parser.add_argument('--num_examples', type=int, required=True)
    parser.add_argument('--num_features', type=int, required=True)
    parser.add_argument('--min_estimators', type=int, default=20)
    parser.add_argument('--max_estimators', type=int, default=200)
    parser.add_argument('--num_estimator_values', type=int, default=3)
    parser.add_argument('--show_plot', action='store_true')
    parser.add_argument('--num_informative', type=int, default=10,
                        help="Used only for make_regression")
    parser.add_argument('--noise', type=float, default=1.0,
                        help="Used only for friedman3")
    parser.add_argument('--dataset', choices=DATASETS.keys(), required=True)
    return parser

if __name__ == '__main__':
    run_simulation(build_parser().parse_args())
