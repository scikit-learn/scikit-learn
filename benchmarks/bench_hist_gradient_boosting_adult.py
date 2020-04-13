import argparse
from time import time

from sklearn.model_selection import train_test_split
from sklearn.datasets import fetch_openml
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.experimental import enable_hist_gradient_boosting  # noqa
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.ensemble._hist_gradient_boosting.utils import (
    get_equivalent_estimator)


parser = argparse.ArgumentParser()
parser.add_argument('--n-leaf-nodes', type=int, default=31)
parser.add_argument('--n-trees', type=int, default=10)
parser.add_argument('--lightgbm', action="store_true", default=False)
parser.add_argument('--learning-rate', type=float, default=1.)
parser.add_argument('--max-bins', type=int, default=255)
parser.add_argument('--no-predict', action="store_true", default=False)
args = parser.parse_args()

n_leaf_nodes = args.n_leaf_nodes
n_trees = args.n_trees
lr = args.learning_rate
max_bins = args.max_bins


def fit(est, data_train, target_train, libname):
    print(f"Fitting a {libname} model...")
    tic = time()
    est.fit(data_train, target_train)
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
    print(f"predicted in {toc - tic:.3f}s, "
          f"ROC AUC: {roc_auc:.4f}, ACC: {acc :.4f}")


data, target = fetch_openml(data_id=179, as_frame=True, return_X_y=True)

# does not support categories in encoding y yet
target = target.cat.codes

n_features = data.shape[1]
n_categorical_features = (data.dtypes == 'category').sum()
n_numerical_features = (data.dtypes == 'float').sum()
print(f"Number of features: {data.shape[1]}")
print(f"Number of categorical featuers: {n_categorical_features}")
print(f"Number of numerical featuers: {n_numerical_features}")

data_train, data_test, target_train, target_test = train_test_split(
    data, target, test_size=.2, random_state=0)


est = HistGradientBoostingClassifier(loss='binary_crossentropy',
                                     learning_rate=lr,
                                     max_iter=n_trees,
                                     max_bins=max_bins,
                                     categorical='pandas',
                                     max_leaf_nodes=n_leaf_nodes,
                                     early_stopping=False,
                                     random_state=0,
                                     verbose=1)
fit(est, data_train, target_train, 'sklearn')
predict(est, data_test, target_test)

# lightgbm infers the cateogires from the dtype
if args.lightgbm:
    est = get_equivalent_estimator(est, lib='lightgbm')
    fit(est, data_train, target_train, 'lightgbm')
    predict(est, data_test, target_test)
