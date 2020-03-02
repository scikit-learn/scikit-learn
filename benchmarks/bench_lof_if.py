"""
=====================
LOF and IF benchmarks
=====================
A test of LocalOutlierFactor (LOF) and IsolationForest (IF) on classical 
anomaly detection datasets. Note that LocalOutlierFactor is not meant to 
predict on a test set and its performance is assessed in an outlier detection 
context:
1. The model is trained on the whole dataset which is assumed to contain
outliers.
2. The ROC curve is computed on the same dataset using the knowledge of the
labels.
In this context there is no need to shuffle the dataset because the model
is trained and tested on the whole dataset. The randomness of this benchmark
is only caused by the random selection of anomalies in the SA dataset.
"""

from time import time
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import LocalOutlierFactor
from sklearn.ensemble import IsolationForest
from sklearn.metrics import roc_curve, auc
from sklearn.datasets import fetch_kddcup99, fetch_covtype, fetch_openml
from sklearn.preprocessing import LabelBinarizer

print(__doc__)

random_state = 1  # to control the random selection of anomalies in SA

# datasets
datasets = ['http', 'smtp', 'SA', 'SF', 'forestcover']
# outlier detection models
models =[('LOF', LocalOutlierFactor(n_neighbors=20, contamination='auto')),
         ('IF', IsolationForest(n_jobs=-1, random_state=random_state,
                                behaviour="new", contamination='auto'))]

plt.figure(figsize=(5,len(datasets)*3))
for dataset_idx, dataset_name in enumerate(datasets):
    plt.subplot(len(datasets), 1, dataset_idx+1)
    # loading and vectorization
    print('loading data: ', str(dataset_idx+1))
    if dataset_name in ['http', 'smtp', 'SA', 'SF']:
        dataset = fetch_kddcup99(subset=dataset_name, percent10=True,
                                 random_state=random_state)
        X = dataset.data
        y = dataset.target

    if dataset_name == 'forestcover':
        dataset = fetch_covtype()
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

    X = X.astype(float)

    print('Estimator processing...')
    for model_name, model in models:
        tstart = time()
        model.fit(X)
        fit_time = time() - tstart
        if model_name == 'LOF':
            scoring = -model.negative_outlier_factor_  # the lower, the more normal
        if model_name == 'IF':
            scoring = -model.fit(X).decision_function(X)

        fpr, tpr, thresholds = roc_curve(y, scoring)
        AUC = auc(fpr, tpr)
        plt.plot(fpr, tpr, lw=1,
                 label=(model_name, ': ROC for %s (area = %0.3f, train-time: %0.2fs)'
                        % (dataset_name, AUC, fit_time)))

    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.legend(loc="lower right")
    if dataset_idx ==0:
        plt.title('Receiver operating characteristic')
    if dataset_idx ==len(datasets)-1:
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')

plt.show()