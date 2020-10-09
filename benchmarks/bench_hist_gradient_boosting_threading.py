from time import time
import argparse
import os
from pprint import pprint

import numpy as np
from threadpoolctl import threadpool_limits
import sklearn
from sklearn.model_selection import train_test_split
# To use this experimental feature, we need to explicitly ask for it:
from sklearn.experimental import enable_hist_gradient_boosting  # noqa
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.datasets import make_classification
from sklearn.datasets import make_regression
from sklearn.ensemble._hist_gradient_boosting.utils import (
    get_equivalent_estimator)


parser = argparse.ArgumentParser()
parser.add_argument('--n-leaf-nodes', type=int, default=31)
parser.add_argument('--n-trees', type=int, default=10)
parser.add_argument('--lightgbm', action="store_true", default=False,
                    help='also benchmark lightgbm')
parser.add_argument('--xgboost', action="store_true", default=False,
                    help='also benchmark xgboost')
parser.add_argument('--catboost', action="store_true", default=False,
                    help='also benchmark catboost')
parser.add_argument('--learning-rate', type=float, default=.1)
parser.add_argument('--problem', type=str, default='classification',
                    choices=['classification', 'regression'])
parser.add_argument('--loss', type=str, default='default')
parser.add_argument('--missing-fraction', type=float, default=0)
parser.add_argument('--n-classes', type=int, default=2)
parser.add_argument('--n-samples', type=int, default=int(1e6))
parser.add_argument('--n-features', type=int, default=100)
parser.add_argument('--max-bins', type=int, default=255)

parser.add_argument('--print-params', action="store_true", default=False)
parser.add_argument('--random-sample-weights', action="store_true",
                    default=False,
                    help="generate and use random sample weights")
parser.add_argument('--plot', action="store_true", default=False,
                    help='show a plot results')
parser.add_argument('--plot-filename', default=None,
                    help='filename to save the figure to disk')
args = parser.parse_args()

n_samples = args.n_samples
n_leaf_nodes = args.n_leaf_nodes
n_trees = args.n_trees
lr = args.learning_rate
max_bins = args.max_bins


print("Data size: %d samples train, %d samples test."
      % (n_samples, n_samples))
print(f"n_features: {args.n_features}")


def get_estimator_and_data():
    if args.problem == 'classification':
        X, y = make_classification(args.n_samples * 2,
                                   n_features=args.n_features,
                                   n_classes=args.n_classes,
                                   n_clusters_per_class=1,
                                   n_informative=args.n_features // 2,
                                   random_state=0)
        return X, y, HistGradientBoostingClassifier
    elif args.problem == 'regression':
        X, y = make_regression(args.n_samples_max * 2,
                               n_features=args.n_features, random_state=0)
        return X, y, HistGradientBoostingRegressor


X, y, Estimator = get_estimator_and_data()
if args.missing_fraction:
    mask = np.random.binomial(1, args.missing_fraction, size=X.shape).astype(
        bool)
    X[mask] = np.nan

if args.random_sample_weights:
    sample_weight = np.random.rand(len(X)) * 10
else:
    sample_weight = None

if sample_weight is not None:
    (X_train_, X_test_, y_train_, y_test_,
     sample_weight_train_, _) = train_test_split(
        X, y, sample_weight, test_size=0.5, random_state=0)
else:
    X_train_, X_test_, y_train_, y_test_ = train_test_split(
        X, y, test_size=0.5, random_state=0)
    sample_weight_train_ = None


sklearn_est = Estimator(
    learning_rate=lr,
    max_iter=n_trees,
    max_bins=max_bins,
    max_leaf_nodes=n_leaf_nodes,
    early_stopping=False,
    random_state=0,
    verbose=0,
)
loss = args.loss
if args.problem == 'classification':
    if loss == 'default':
        # loss='auto' does not work with get_equivalent_estimator()
        loss = 'binary_crossentropy' if args.n_classes == 2 else \
            'categorical_crossentropy'
else:
    # regression
    if loss == 'default':
        loss = 'least_squares'
sklearn_est.set_params(loss=loss)


if args.print_params:
    print("scikit-learn")
    pprint(sklearn_est.get_params())

    for libname in ["lightgbm", "xgboost", "catboost"]:
        if getattr(args, libname):
            print(libname)
            est = get_equivalent_estimator(sklearn_est, lib=libname)
            pprint(est.get_params())


def one_run(n_threads, n_samples):
    X_train = X_train_[:n_samples]
    X_test = X_test_[:n_samples]
    y_train = y_train_[:n_samples]
    y_test = y_test_[:n_samples]
    if sample_weight is not None:
        sample_weight_train = sample_weight_train_[:n_samples]
    else:
        sample_weight_train = None
    assert X_train.shape[0] == n_samples
    assert X_test.shape[0] == n_samples
    print("Fitting a sklearn model...")
    tic = time()
    est = sklearn.base.clone(sklearn_est)

    with threadpool_limits(n_threads, user_api="openmp"):
        est.fit(X_train, y_train, sample_weight=sample_weight_train)
        sklearn_fit_duration = time() - tic
        tic = time()
        sklearn_score = est.score(X_test, y_test)
        sklearn_score_duration = time() - tic
    print("score: {:.4f}".format(sklearn_score))
    print("fit duration: {:.3f}s,".format(sklearn_fit_duration))
    print("score duration: {:.3f}s,".format(sklearn_score_duration))

    lightgbm_score = None
    lightgbm_fit_duration = None
    lightgbm_score_duration = None
    if args.lightgbm:
        print("Fitting a LightGBM model...")
        lightgbm_est = get_equivalent_estimator(est, lib='lightgbm')
        lightgbm_est.set_params(num_threads=n_threads)

        tic = time()
        lightgbm_est.fit(X_train, y_train, sample_weight=sample_weight_train)
        lightgbm_fit_duration = time() - tic
        tic = time()
        lightgbm_score = lightgbm_est.score(X_test, y_test)
        lightgbm_score_duration = time() - tic
        print("score: {:.4f}".format(lightgbm_score))
        print("fit duration: {:.3f}s,".format(lightgbm_fit_duration))
        print("score duration: {:.3f}s,".format(lightgbm_score_duration))

    xgb_score = None
    xgb_fit_duration = None
    xgb_score_duration = None
    if args.xgboost:
        print("Fitting an XGBoost model...")
        xgb_est = get_equivalent_estimator(est, lib='xgboost')
        xgb_est.set_params(nthread=n_threads)

        tic = time()
        xgb_est.fit(X_train, y_train, sample_weight=sample_weight_train)
        xgb_fit_duration = time() - tic
        tic = time()
        xgb_score = xgb_est.score(X_test, y_test)
        xgb_score_duration = time() - tic
        print("score: {:.4f}".format(xgb_score))
        print("fit duration: {:.3f}s,".format(xgb_fit_duration))
        print("score duration: {:.3f}s,".format(xgb_score_duration))

    cat_score = None
    cat_fit_duration = None
    cat_score_duration = None
    if args.catboost:
        print("Fitting a CatBoost model...")
        cat_est = get_equivalent_estimator(est, lib='catboost')
        cat_est.set_params(thread_count=n_threads)

        tic = time()
        cat_est.fit(X_train, y_train, sample_weight=sample_weight_train)
        cat_fit_duration = time() - tic
        tic = time()
        cat_score = cat_est.score(X_test, y_test)
        cat_score_duration = time() - tic
        print("score: {:.4f}".format(cat_score))
        print("fit duration: {:.3f}s,".format(cat_fit_duration))
        print("score duration: {:.3f}s,".format(cat_score_duration))

    return (sklearn_score, sklearn_fit_duration, sklearn_score_duration,
            lightgbm_score, lightgbm_fit_duration, lightgbm_score_duration,
            xgb_score, xgb_fit_duration, xgb_score_duration,
            cat_score, cat_fit_duration, cat_score_duration)


max_threads = os.cpu_count()
n_threads_list = [2 ** i for i in range(8) if (2 ** i) < max_threads]
n_threads_list.append(max_threads)

sklearn_scores = []
sklearn_fit_durations = []
sklearn_score_durations = []
lightgbm_scores = []
lightgbm_fit_durations = []
lightgbm_score_durations = []
xgb_scores = []
xgb_fit_durations = []
xgb_score_durations = []
cat_scores = []
cat_fit_durations = []
cat_score_durations = []

for n_threads in n_threads_list:
    print(f"n_threads: {n_threads}")
    (
        sklearn_score,
        sklearn_fit_duration,
        sklearn_score_duration,
        lightgbm_score,
        lightgbm_fit_duration,
        lightgbm_score_duration,
        xgb_score,
        xgb_fit_duration,
        xgb_score_duration,
        cat_score,
        cat_fit_duration,
        cat_score_duration
    ) = one_run(n_threads, n_samples)

    for scores, score in (
            (sklearn_scores, sklearn_score),
            (sklearn_fit_durations, sklearn_fit_duration),
            (sklearn_score_durations, sklearn_score_duration),
            (lightgbm_scores, lightgbm_score),
            (lightgbm_fit_durations, lightgbm_fit_duration),
            (lightgbm_score_durations, lightgbm_score_duration),
            (xgb_scores, xgb_score),
            (xgb_fit_durations, xgb_fit_duration),
            (xgb_score_durations, xgb_score_duration),
            (cat_scores, cat_score),
            (cat_fit_durations, cat_fit_duration),
            (cat_score_durations, cat_score_duration)):
        scores.append(score)


if args.plot or args.plot_filename:
    import matplotlib.pyplot as plt
    import matplotlib

    fig, axs = plt.subplots(2, figsize=(12, 12))

    label = f"sklearn {sklearn.__version__}"
    axs[0].plot(n_threads_list, sklearn_fit_durations, label=label)
    axs[1].plot(n_threads_list, sklearn_score_durations, label=label)

    if args.lightgbm:
        import lightgbm
        label = f'LightGBM {lightgbm.__version__}'
        axs[0].plot(n_threads_list, lightgbm_fit_durations, label=label)
        axs[1].plot(n_threads_list, lightgbm_score_durations, label=label)

    if args.xgboost:
        import xgboost
        label = f'XGBoost {xgboost.__version__}'
        axs[0].plot(n_threads_list, xgb_fit_durations, label=label)
        axs[1].plot(n_threads_list, xgb_score_durations, label=label)

    if args.catboost:
        import catboost
        label = f'CatBoost {catboost.__version__}'
        axs[0].plot(n_threads_list, cat_fit_durations, label=label)
        axs[1].plot(n_threads_list, cat_score_durations, label=label)

    for ax in axs:
        ax.set_xscale('log')
        ax.set_xlabel('n_threads')
        ax.set_ylabel('duration (s)')
        ax.set_ylim(0, None)
        ax.set_xticks(n_threads_list)
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.legend(loc='best')

    axs[0].set_title('fit duration (s)')
    axs[1].set_title('score duration (s)')

    title = args.problem
    if args.problem == 'classification':
        title += ' n_classes = {}'.format(args.n_classes)
    fig.suptitle(title)

    plt.tight_layout()

    if args.plot_filename:
        plt.savefig(args.plot_filename)

    if args.plot:
        plt.show()
