import argparse
from time import time

import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.compose import make_column_transformer, make_column_selector
from sklearn.datasets import fetch_openml
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.ensemble._hist_gradient_boosting.utils import get_equivalent_estimator
from sklearn.preprocessing import OrdinalEncoder


parser = argparse.ArgumentParser()
parser.add_argument("--n-leaf-nodes", type=int, default=31)
parser.add_argument("--n-trees", type=int, default=100)
parser.add_argument("--lightgbm", action="store_true", default=False)
parser.add_argument("--learning-rate", type=float, default=0.1)
parser.add_argument("--max-bins", type=int, default=255)
parser.add_argument("--no-predict", action="store_true", default=False)
parser.add_argument("--verbose", action="store_true", default=False)
args = parser.parse_args()

n_leaf_nodes = args.n_leaf_nodes
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


def predict(est, data_test, target_test):
    if args.no_predict:
        return
    tic = time()
    predicted_test = est.predict(data_test)
    predicted_proba_test = est.predict_proba(data_test)
    toc = time()
    roc_auc = roc_auc_score(target_test, predicted_proba_test[:, 1])
    acc = accuracy_score(target_test, predicted_test)
    print(f"predicted in {toc - tic:.3f}s, ROC AUC: {roc_auc:.4f}, ACC: {acc :.4f}")


data = fetch_openml(data_id=179, as_frame=True, parser="pandas")  # adult dataset
X, y = data.data, data.target

# Ordinal encode the categories to use the native support available in HGBDT
cat_columns = make_column_selector(dtype_include="category")(X)
preprocessing = make_column_transformer(
    (OrdinalEncoder(), cat_columns),
    remainder="passthrough",
    verbose_feature_names_out=False,
)
X = pd.DataFrame(
    preprocessing.fit_transform(X),
    columns=preprocessing.get_feature_names_out(),
)

n_classes = len(np.unique(y))
n_features = X.shape[1]
n_categorical_features = len(cat_columns)
n_numerical_features = n_features - n_categorical_features
print(f"Number of features: {n_features}")
print(f"Number of categorical features: {n_categorical_features}")
print(f"Number of numerical features: {n_numerical_features}")

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)

is_categorical = [True] * n_categorical_features + [False] * n_numerical_features
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

fit(est, X_train, y_train, "sklearn")
predict(est, X_test, y_test)

if args.lightgbm:
    est = get_equivalent_estimator(est, lib="lightgbm", n_classes=n_classes)
    est.set_params(max_cat_to_onehot=1)  # dont use OHE
    categorical_features = [
        f_idx for (f_idx, is_cat) in enumerate(is_categorical) if is_cat
    ]
    fit(est, X_train, y_train, "lightgbm", categorical_feature=categorical_features)
    predict(est, X_test, y_test)
