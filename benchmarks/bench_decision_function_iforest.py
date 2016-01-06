"""
=============================================
IsolationForest's decision function benchmark
=============================================

An histogram of IsolationForest's decision function on classical anomaly
detection datasets.

This is an example which plots the decision function of IsolationForest
for normal data and for outliers on real-world datasets.

"""
print(__doc__)

import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import IsolationForest
from sklearn.datasets import fetch_kddcup99, fetch_covtype, fetch_mldata
from sklearn.preprocessing import LabelBinarizer
from sklearn.utils import shuffle as sh

np.random.seed(1)


datasets = ['http', 'smtp', 'SA', 'SF', 'shuttle', 'forestcover']

for dat_nb, dat in enumerate(datasets):
    # loading and vectorization
    print('loading data')
    if dat in ['http', 'smtp', 'SA', 'SF']:
        dataset = fetch_kddcup99(subset=dat, shuffle=True, percent10=True)
        X = dataset.data
        y = dataset.target

    if dat == 'shuttle':
        dataset = fetch_mldata('shuttle')
        X = dataset.data
        y = dataset.target
        sh(X, y)
        # we remove data with label 4
        # normal data are then those of class 1
        s = (y != 4)
        X = X[s, :]
        y = y[s]
        y = (y != 1).astype(int)

    if dat == 'forestcover':
        dataset = fetch_covtype(shuffle=True)
        X = dataset.data
        y = dataset.target
        # normal data are those with attribute 2
        # abnormal those with attribute 4
        s = (y == 2) + (y == 4)
        X = X[s, :]
        y = y[s]
        y = (y != 2).astype(int)

    print('vectorizing data')

    if dat == 'SF':
        lb = LabelBinarizer()
        lb.fit(X[:, 1])
        x1 = lb.transform(X[:, 1])
        X = np.c_[X[:, :1], x1, X[:, 2:]]
        y = (y != 'normal.').astype(int)

    if dat == 'SA':
        lb = LabelBinarizer()
        lb.fit(X[:, 1])
        x1 = lb.transform(X[:, 1])
        lb.fit(X[:, 2])
        x2 = lb.transform(X[:, 2])
        lb.fit(X[:, 3])
        x3 = lb.transform(X[:, 3])
        X = np.c_[X[:, :1], x1, x2, x3, X[:, 4:]]
        y = (y != 'normal.').astype(int)

    if dat == 'http' or dat == 'smtp':
        y = (y != 'normal.').astype(int)

    n_samples, n_features = np.shape(X)
    n_samples_train = n_samples // 2
    n_samples_test = n_samples - n_samples_train

    X = X.astype(float)
    X_train = X[:n_samples_train, :]
    X_test = X[n_samples_train:, :]
    y_train = y[:n_samples_train]
    y_test = y[n_samples_train:]

    print('IsolationForest processing...')
    model = IsolationForest(n_jobs=-1)
    model.fit(X_train)

    scoring = 0.5 + model.decision_function(X_test)  # to center around 0

    f, ax = plt.subplots(3, sharex=True, sharey=True)
    ax[0].hist(scoring, np.linspace(-0.5, 0.5, 200), color='black')
    ax[0].set_title('decision function for %s dataset' % dat, size=20)
    ax[0].legend(loc="lower right")
    ax[1].hist(scoring[y_test == 0], np.linspace(-0.5, 0.5, 200), color='b',
               label='normal data')
    ax[1].legend(loc="lower right")
    ax[2].hist(scoring[y_test == 1], np.linspace(-0.5, 0.5, 200), color='r',
               label='outliers')
    ax[2].legend(loc="lower right")
    plt.show()
