"""
=====================================
SGDOneClassSVM benchmark
=====================================
This benchmark compares the :class:`SGDOneClassSVM` with :class:`OneClassSVM`.
The former is an online One-Class SVM implemented with a Stochastic Gradient
Descent (SGD). The latter is based on the LibSVM implementation. The
complexity of :class:`SGDOneClassSVM` is linear in the number of samples
whereas the one of :class:`OneClassSVM` is at best quadratic in the number of
samples. We here compare the performance in terms of AUC and training time on
classical anomaly detection datasets.

The :class:`OneClassSVM` is applied with a Gaussian kernel and we therefore
use a kernel approximation prior to the application of :class:`SGDOneClassSVM`.
"""

from time import time
import numpy as np

from scipy.interpolate import interp1d

from sklearn.metrics import roc_curve, auc
from sklearn.datasets import fetch_kddcup99, fetch_covtype
from sklearn.preprocessing import LabelBinarizer, StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.utils import shuffle
from sklearn.kernel_approximation import Nystroem
from sklearn.svm import OneClassSVM
from sklearn.linear_model import SGDOneClassSVM

import matplotlib.pyplot as plt
import matplotlib

font = {"weight": "normal", "size": 15}

matplotlib.rc("font", **font)

print(__doc__)


def print_outlier_ratio(y):
    """
    Helper function to show the distinct value count of element in the target.
    Useful indicator for the datasets used in bench_isolation_forest.py.
    """
    uniq, cnt = np.unique(y, return_counts=True)
    print("----- Target count values: ")
    for u, c in zip(uniq, cnt):
        print("------ %s -> %d occurrences" % (str(u), c))
    print("----- Outlier ratio: %.5f" % (np.min(cnt) / len(y)))


# for roc curve computation
n_axis = 1000
x_axis = np.linspace(0, 1, n_axis)

datasets = ["http", "smtp", "SA", "SF", "forestcover"]

novelty_detection = False  # if False, training set polluted by outliers

random_states = [42]
nu = 0.05

results_libsvm = np.empty((len(datasets), n_axis + 5))
results_online = np.empty((len(datasets), n_axis + 5))

for dat, dataset_name in enumerate(datasets):

    print(dataset_name)

    # Loading datasets
    if dataset_name in ["http", "smtp", "SA", "SF"]:
        dataset = fetch_kddcup99(
            subset=dataset_name, shuffle=False, percent10=False, random_state=88
        )
        X = dataset.data
        y = dataset.target

    if dataset_name == "forestcover":
        dataset = fetch_covtype(shuffle=False)
        X = dataset.data
        y = dataset.target
        # normal data are those with attribute 2
        # abnormal those with attribute 4
        s = (y == 2) + (y == 4)
        X = X[s, :]
        y = y[s]
        y = (y != 2).astype(int)

    # Vectorizing data
    if dataset_name == "SF":
        # Casting type of X (object) as string is needed for string categorical
        # features to apply LabelBinarizer
        lb = LabelBinarizer()
        x1 = lb.fit_transform(X[:, 1].astype(str))
        X = np.c_[X[:, :1], x1, X[:, 2:]]
        y = (y != b"normal.").astype(int)

    if dataset_name == "SA":
        lb = LabelBinarizer()
        # Casting type of X (object) as string is needed for string categorical
        # features to apply LabelBinarizer
        x1 = lb.fit_transform(X[:, 1].astype(str))
        x2 = lb.fit_transform(X[:, 2].astype(str))
        x3 = lb.fit_transform(X[:, 3].astype(str))
        X = np.c_[X[:, :1], x1, x2, x3, X[:, 4:]]
        y = (y != b"normal.").astype(int)

    if dataset_name in ["http", "smtp"]:
        y = (y != b"normal.").astype(int)

    print_outlier_ratio(y)

    n_samples, n_features = np.shape(X)
    if dataset_name == "SA":  # LibSVM too long with n_samples // 2
        n_samples_train = n_samples // 20
    else:
        n_samples_train = n_samples // 2

    n_samples_test = n_samples - n_samples_train
    print("n_train: ", n_samples_train)
    print("n_features: ", n_features)

    tpr_libsvm = np.zeros(n_axis)
    tpr_online = np.zeros(n_axis)
    fit_time_libsvm = 0
    fit_time_online = 0
    predict_time_libsvm = 0
    predict_time_online = 0

    X = X.astype(float)

    gamma = 1 / n_features  # OCSVM default parameter

    for random_state in random_states:

        print("random state: %s" % random_state)

        X, y = shuffle(X, y, random_state=random_state)
        X_train = X[:n_samples_train]
        X_test = X[n_samples_train:]
        y_train = y[:n_samples_train]
        y_test = y[n_samples_train:]

        if novelty_detection:
            X_train = X_train[y_train == 0]
            y_train = y_train[y_train == 0]

        std = StandardScaler()

        print("----------- LibSVM OCSVM ------------")
        ocsvm = OneClassSVM(kernel="rbf", gamma=gamma, nu=nu)
        pipe_libsvm = make_pipeline(std, ocsvm)

        tstart = time()
        pipe_libsvm.fit(X_train)
        fit_time_libsvm += time() - tstart

        tstart = time()
        # scoring such that the lower, the more normal
        scoring = -pipe_libsvm.decision_function(X_test)
        predict_time_libsvm += time() - tstart
        fpr_libsvm_, tpr_libsvm_, _ = roc_curve(y_test, scoring)

        f_libsvm = interp1d(fpr_libsvm_, tpr_libsvm_)
        tpr_libsvm += f_libsvm(x_axis)

        print("----------- Online OCSVM ------------")
        nystroem = Nystroem(gamma=gamma, random_state=random_state)
        online_ocsvm = SGDOneClassSVM(nu=nu, random_state=random_state)
        pipe_online = make_pipeline(std, nystroem, online_ocsvm)

        tstart = time()
        pipe_online.fit(X_train)
        fit_time_online += time() - tstart

        tstart = time()
        # scoring such that the lower, the more normal
        scoring = -pipe_online.decision_function(X_test)
        predict_time_online += time() - tstart
        fpr_online_, tpr_online_, _ = roc_curve(y_test, scoring)

        f_online = interp1d(fpr_online_, tpr_online_)
        tpr_online += f_online(x_axis)

    tpr_libsvm /= len(random_states)
    tpr_libsvm[0] = 0.0
    fit_time_libsvm /= len(random_states)
    predict_time_libsvm /= len(random_states)
    auc_libsvm = auc(x_axis, tpr_libsvm)

    results_libsvm[dat] = [
        fit_time_libsvm,
        predict_time_libsvm,
        auc_libsvm,
        n_samples_train,
        n_features,
    ] + list(tpr_libsvm)

    tpr_online /= len(random_states)
    tpr_online[0] = 0.0
    fit_time_online /= len(random_states)
    predict_time_online /= len(random_states)
    auc_online = auc(x_axis, tpr_online)

    results_online[dat] = [
        fit_time_online,
        predict_time_online,
        auc_online,
        n_samples_train,
        n_features,
    ] + list(tpr_libsvm)


# -------- Plotting bar charts -------------
fit_time_libsvm_all = results_libsvm[:, 0]
predict_time_libsvm_all = results_libsvm[:, 1]
auc_libsvm_all = results_libsvm[:, 2]
n_train_all = results_libsvm[:, 3]
n_features_all = results_libsvm[:, 4]

fit_time_online_all = results_online[:, 0]
predict_time_online_all = results_online[:, 1]
auc_online_all = results_online[:, 2]


width = 0.7
ind = 2 * np.arange(len(datasets))
x_tickslabels = [
    (name + "\n" + r"$n={:,d}$" + "\n" + r"$d={:d}$").format(int(n), int(d))
    for name, n, d in zip(datasets, n_train_all, n_features_all)
]


def autolabel_auc(rects, ax):
    """Attach a text label above each bar displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.text(
            rect.get_x() + rect.get_width() / 2.0,
            1.05 * height,
            "%.3f" % height,
            ha="center",
            va="bottom",
        )


def autolabel_time(rects, ax):
    """Attach a text label above each bar displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.text(
            rect.get_x() + rect.get_width() / 2.0,
            1.05 * height,
            "%.1f" % height,
            ha="center",
            va="bottom",
        )


fig, ax = plt.subplots(figsize=(15, 8))
ax.set_ylabel("AUC")
ax.set_ylim((0, 1.3))
rect_libsvm = ax.bar(ind, auc_libsvm_all, width=width, color="r")
rect_online = ax.bar(ind + width, auc_online_all, width=width, color="y")
ax.legend((rect_libsvm[0], rect_online[0]), ("LibSVM", "Online SVM"))
ax.set_xticks(ind + width / 2)
ax.set_xticklabels(x_tickslabels)
autolabel_auc(rect_libsvm, ax)
autolabel_auc(rect_online, ax)
plt.show()


fig, ax = plt.subplots(figsize=(15, 8))
ax.set_ylabel("Training time (sec) - Log scale")
ax.set_yscale("log")
rect_libsvm = ax.bar(ind, fit_time_libsvm_all, color="r", width=width)
rect_online = ax.bar(ind + width, fit_time_online_all, color="y", width=width)
ax.legend((rect_libsvm[0], rect_online[0]), ("LibSVM", "Online SVM"))
ax.set_xticks(ind + width / 2)
ax.set_xticklabels(x_tickslabels)
autolabel_time(rect_libsvm, ax)
autolabel_time(rect_online, ax)
plt.show()


fig, ax = plt.subplots(figsize=(15, 8))
ax.set_ylabel("Testing time (sec) - Log scale")
ax.set_yscale("log")
rect_libsvm = ax.bar(ind, predict_time_libsvm_all, color="r", width=width)
rect_online = ax.bar(ind + width, predict_time_online_all, color="y", width=width)
ax.legend((rect_libsvm[0], rect_online[0]), ("LibSVM", "Online SVM"))
ax.set_xticks(ind + width / 2)
ax.set_xticklabels(x_tickslabels)
autolabel_time(rect_libsvm, ax)
autolabel_time(rect_online, ax)
plt.show()
