from urllib.request import urlretrieve
import os
from gzip import GzipFile
from time import time
import argparse

import numpy as np
import pandas as pd
from joblib import Memory
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.ensemble._hist_gradient_boosting.utils import get_equivalent_estimator


parser = argparse.ArgumentParser()
parser.add_argument("--n-leaf-nodes", type=int, default=31)
parser.add_argument("--n-trees", type=int, default=10)
parser.add_argument("--lightgbm", action="store_true", default=False)
parser.add_argument("--xgboost", action="store_true", default=False)
parser.add_argument("--catboost", action="store_true", default=False)
parser.add_argument("--learning-rate", type=float, default=1.0)
parser.add_argument("--subsample", type=int, default=None)
parser.add_argument("--max-bins", type=int, default=255)
parser.add_argument("--no-predict", action="store_true", default=False)
parser.add_argument("--cache-loc", type=str, default="/tmp")
args = parser.parse_args()

HERE = os.path.dirname(__file__)
URL = "https://archive.ics.uci.edu/ml/machine-learning-databases/00280/HIGGS.csv.gz"
m = Memory(location=args.cache_loc, mmap_mode="r")

n_leaf_nodes = args.n_leaf_nodes
n_trees = args.n_trees
subsample = args.subsample
lr = args.learning_rate
max_bins = args.max_bins


@m.cache
def load_data():
    filename = os.path.join(HERE, URL.rsplit("/", 1)[-1])
    if not os.path.exists(filename):
        print(f"Downloading {URL} to {filename} (2.6 GB)...")
        urlretrieve(URL, filename)
        print("done.")

    print(f"Parsing {filename}...")
    tic = time()
    with GzipFile(filename) as f:
        df = pd.read_csv(f, header=None, dtype=np.float32)
    toc = time()
    print(f"Loaded {df.values.nbytes / 1e9:0.3f} GB in {toc - tic:0.3f}s")
    return df


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
    print(f"predicted in {toc - tic:.3f}s, ROC AUC: {roc_auc:.4f}, ACC: {acc :.4f}")


df = load_data()
target = df.values[:, 0]
data = np.ascontiguousarray(df.values[:, 1:])
data_train, data_test, target_train, target_test = train_test_split(
    data, target, test_size=0.2, random_state=0
)
n_classes = len(np.unique(target))

if subsample is not None:
    data_train, target_train = data_train[:subsample], target_train[:subsample]

n_samples, n_features = data_train.shape
print(f"Training set with {n_samples} records with {n_features} features.")

est = HistGradientBoostingClassifier(
    loss="log_loss",
    learning_rate=lr,
    max_iter=n_trees,
    max_bins=max_bins,
    max_leaf_nodes=n_leaf_nodes,
    early_stopping=False,
    random_state=0,
    verbose=1,
)
fit(est, data_train, target_train, "sklearn")
predict(est, data_test, target_test)

if args.lightgbm:
    est = get_equivalent_estimator(est, lib="lightgbm", n_classes=n_classes)
    fit(est, data_train, target_train, "lightgbm")
    predict(est, data_test, target_test)

if args.xgboost:
    est = get_equivalent_estimator(est, lib="xgboost", n_classes=n_classes)
    fit(est, data_train, target_train, "xgboost")
    predict(est, data_test, target_test)

if args.catboost:
    est = get_equivalent_estimator(est, lib="catboost", n_classes=n_classes)
    fit(est, data_train, target_train, "catboost")
    predict(est, data_test, target_test)
