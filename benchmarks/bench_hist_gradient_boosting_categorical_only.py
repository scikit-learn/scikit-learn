import argparse
from time import time

from sklearn.datasets import make_classification
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.ensemble._hist_gradient_boosting.utils import get_equivalent_estimator
from sklearn.preprocessing import KBinsDiscretizer

parser = argparse.ArgumentParser()
parser.add_argument("--n-leaf-nodes", type=int, default=31)
parser.add_argument("--n-trees", type=int, default=100)
parser.add_argument("--n-features", type=int, default=20)
parser.add_argument("--n-cats", type=int, default=20)
parser.add_argument("--n-samples", type=int, default=10_000)
parser.add_argument("--lightgbm", action="store_true", default=False)
parser.add_argument("--learning-rate", type=float, default=0.1)
parser.add_argument("--max-bins", type=int, default=255)
parser.add_argument("--no-predict", action="store_true", default=False)
parser.add_argument("--verbose", action="store_true", default=False)
args = parser.parse_args()

n_leaf_nodes = args.n_leaf_nodes
n_features = args.n_features
n_categories = args.n_cats
n_samples = args.n_samples
n_trees = args.n_trees
lr = args.learning_rate
max_bins = args.max_bins
verbose = args.verbose


def fit(est, data_train, target_train, libname, **fit_params):
    print(f"Fitting a {libname} model...")
    tic = time()
    est.fit(data_train, target_train, **fit_params)
    toc = time()
    print(f"fitted in {toc - tic:.3f}s")


def predict(est, data_test):
    # We don't report accuracy or ROC because the dataset doesn't really make
    # sense: we treat ordered features as un-ordered categories.
    if args.no_predict:
        return
    tic = time()
    est.predict(data_test)
    toc = time()
    print(f"predicted in {toc - tic:.3f}s")


X, y = make_classification(n_samples=n_samples, n_features=n_features, random_state=0)

X = KBinsDiscretizer(n_bins=n_categories, encode="ordinal").fit_transform(X)

print(f"Number of features: {n_features}")
print(f"Number of samples: {n_samples}")

is_categorical = [True] * n_features
est = HistGradientBoostingClassifier(
    loss="log_loss",
    learning_rate=lr,
    max_iter=n_trees,
    max_bins=max_bins,
    max_leaf_nodes=n_leaf_nodes,
    categorical_features=is_categorical,
    early_stopping=False,
    random_state=0,
    verbose=verbose,
)

fit(est, X, y, "sklearn")
predict(est, X)

if args.lightgbm:
    est = get_equivalent_estimator(est, lib="lightgbm", n_classes=2)
    est.set_params(max_cat_to_onehot=1)  # dont use OHE
    categorical_features = list(range(n_features))
    fit(est, X, y, "lightgbm", categorical_feature=categorical_features)
    predict(est, X)
