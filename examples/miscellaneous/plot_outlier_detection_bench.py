"""
==========================================
Evaluation of outlier detection estimators
==========================================

This example benchmarks outlier detection algorithms, Local Outlier
Factor (LOF) and Isolation Forest (IForest), using ROC curves on
classical anomaly detection datasets. The algorithm performance
is assessed in an outlier detection context:

1. The algorithms are trained on the whole dataset which is assumed to
contain outliers.

2. The ROC curve is computed on the same dataset using the knowledge
of the labels.

There is no particular reason to choose algorithms LOF and ROC. The goal
is to show that different algorithm performs well on different datasets.

Interpreting the ROC plot
-------------------------
The algorithm performance relates to how good the true positive rate (TPR)
is at low value of the false positive rate (FPR). The better algorithms
have the curve on the top-left of the plot and the area under curve (AUC)
close to 1. The diagonal dashed line represents a random classification
of outliers and inliers.
"""

from time import time
import math
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import LocalOutlierFactor
from sklearn.ensemble import IsolationForest
from sklearn.metrics import roc_curve, auc
from sklearn.datasets import fetch_kddcup99, fetch_covtype, fetch_openml
from sklearn.preprocessing import LabelBinarizer

print(__doc__)

rng = np.random.RandomState(42)

# datasets
datasets = [
	"http",
	"smtp",
	"SA",
	"SF",
	"forestcover",
	"glass",
	"wdbc",
	"cardiotocography"
]

# outlier detection models
models = [
    ("LOF", LocalOutlierFactor(n_neighbors=20, contamination="auto")),
    ("IForest", IsolationForest(random_state=rng, contamination="auto")),
]

rows = math.ceil(len(datasets) / 2)
plt.figure(figsize=(10, rows * 3))
for dataset_idx, dataset_name in enumerate(datasets):
    plt.subplot(rows, 2, dataset_idx + 1)

    # loading and vectorization
    print(f"Loading {dataset_name} data")
    if dataset_name in ["http", "smtp", "SA", "SF"]:
        dataset = fetch_kddcup99(subset=dataset_name, percent10=True, random_state=rng)
        X = dataset.data
        y = dataset.target

    print("vectorizing data")
    if dataset_name == "SF":
        idx = rng.choice(X.shape[0], int(X.shape[0] * 0.1), replace=False)
        X = X[idx]  # reduce the sample size to speed up computation
        y = y[idx]
        lb = LabelBinarizer()
        x1 = lb.fit_transform(X[:, 1].astype(str))
        X = np.c_[X[:, :1], x1, X[:, 2:]]
        y = (y != b"normal.").astype(int)

    if dataset_name == "SA":
        idx = rng.choice(X.shape[0], int(X.shape[0] * 0.1), replace=False)
        X = X[idx]  # reduce the sample size to speed up computation
        y = y[idx]
        lb = LabelBinarizer()
        x1 = lb.fit_transform(X[:, 1].astype(str))
        x2 = lb.fit_transform(X[:, 2].astype(str))
        x3 = lb.fit_transform(X[:, 3].astype(str))
        X = np.c_[X[:, :1], x1, x2, x3, X[:, 4:]]
        y = (y != b"normal.").astype(int)

    if dataset_name == "http" or dataset_name == "smtp":
        y = (y != b"normal.").astype(int)

    if dataset_name == "forestcover":
        dataset = fetch_covtype()
        X = dataset.data
        y = dataset.target
        idx = rng.choice(X.shape[0], int(X.shape[0] * 0.1), replace=False)
        X = X[idx]  # reduce the sample size to speed up computation
        y = y[idx]

        # normal data are those with attribute 2
        # abnormal those with attribute 4
        s = (y == 2) + (y == 4)
        X = X[s, :]
        y = y[s]
        y = (y != 2).astype(int)

    if dataset_name in ["glass", "wdbc", "cardiotocography"]:
        dataset = fetch_openml(name=dataset_name, version=1,as_frame=False)
        X = dataset.data
        y = dataset.target

        if dataset_name == "glass":
            s = y == "tableware"
            y = s.astype(int)

        if dataset_name == "wdbc":
            s = y == "2"
            y = s.astype(int)
            X_mal, y_mal = X[s], y[s]
            X_ben, y_ben = X[~s], y[~s]

            # downsampled to 39 points (9.8% outliers) to speed up computation
            idx = rng.choice(y_mal.shape[0], 39, replace=False)
            X_mal2 = X_mal[idx]
            y_mal2 = y_mal[idx]
            X = np.concatenate((X_ben, X_mal2), axis=0)
            y = np.concatenate((y_ben, y_mal2), axis=0)

        if dataset_name == "cardiotocography":
            s = y == "3"
            y = s.astype(int)

    print("Estimator processing...")
    for model_name, model in models:
        tstart = time()
        model.fit(X)
        fit_time = time() - tstart

        # use oppositve of default scoring so that the lower score, the more
        # normal
        if model_name == "LOF":
            scoring = -model.negative_outlier_factor_

        if model_name == "IForest":
            scoring = -model.fit(X).decision_function(X)

        fpr, tpr, thresholds = roc_curve(y, scoring)
        area = auc(fpr, tpr)
        label_ = (f"{model_name} (AUC = {area:0.3f}, train time: {fit_time:0.2f})")
        plt.plot(fpr, tpr, lw=1, label=label_)

    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.legend(loc="lower right")
    plt.title(dataset_name)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
plt.tight_layout(pad=2.0)  # spacing between subplots
plt.show()
