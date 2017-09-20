"""
============================
LocalOutlierFactor benchmark
============================

A test of LocalOutlierFactor on classical anomaly detection datasets.

Note that LocalOutlierFactor is not meant to predict on a test set and its
performance is assessed in an outlier detection context:
1. The model is trained on a dataset containing outliers.
2. The ROC curve is computed on the whole dataset using the knowledge of the
labels.

"""

from time import time
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import LocalOutlierFactor
from sklearn.metrics import roc_curve, auc
from sklearn.datasets import fetch_kddcup99, fetch_covtype, fetch_mldata
from sklearn.preprocessing import LabelBinarizer
from sklearn.utils import shuffle as sh

print(__doc__)

np.random.seed(2)

# datasets available: ['http', 'smtp', 'SA', 'SF', 'shuttle', 'forestcover']
datasets = ['http', 'smtp', 'SA', 'SF', 'shuttle', 'forestcover']

for dataset_name in datasets:
    # loading and vectorization
    print('loading data')
    if dataset_name in ['http', 'smtp', 'SA', 'SF']:
        dataset = fetch_kddcup99(subset=dataset_name, shuffle=True,
                                 percent10=False)
        X = dataset.data
        y = dataset.target

    if dataset_name == 'shuttle':
        dataset = fetch_mldata('shuttle')
        X = dataset.data
        y = dataset.target
        X, y = sh(X, y)
        # we remove data with label 4
        # normal data are then those of class 1
        s = (y != 4)
        X = X[s, :]
        y = y[s]
        y = (y != 1).astype(int)

    if dataset_name == 'forestcover':
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

    if dataset_name == 'SF':
        lb = LabelBinarizer()
        x1 = lb.fit_transform(X[:, 1].astype(str))
        X = np.c_[X[:, :1], x1, X[:, 2:]]
        y = (y != b'normal.').astype(int)

    if dataset_name == 'SA':
        lb = LabelBinarizer()
        x1 = lb.fit_transform(X[:, 1].astype(str))
        x2 = lb.fit_transform(X[:, 2].astype(str))
        x3 = lb.fit_transform(X[:, 3].astype(str))
        X = np.c_[X[:, :1], x1, x2, x3, X[:, 4:]]
        y = (y != b'normal.').astype(int)

    if dataset_name == 'http' or dataset_name == 'smtp':
        y = (y != b'normal.').astype(int)

    print('LocalOutlierFactor processing...')
    model = LocalOutlierFactor(n_neighbors=20)
    tstart = time()
    model.fit(X)
    fit_time = time() - tstart
    tstart = time()

    scoring = -model.negative_outlier_factor_  # the lower, the more normal
    predict_time = time() - tstart
    fpr, tpr, thresholds = roc_curve(y, scoring)
    AUC = auc(fpr, tpr)
    plt.plot(fpr, tpr, lw=1,
             label=('ROC for %s (area = %0.3f, train-time: %0.2fs,'
                    'test-time: %0.2fs)' % (dataset_name, AUC, fit_time,
                                            predict_time)))

plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic')
plt.legend(loc="lower right")
plt.show()
