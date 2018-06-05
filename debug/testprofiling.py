import numpy as np
from scipy.stats import bernoulli
from sklearn.metrics.classification import precision_recall_fscore_support, \
    precision_recall_fscore_support_with_multilabel_confusion_matrix

n_samples = 30000
n_labels = 2000
y_true = bernoulli.rvs(np.ones((n_samples, n_labels)) / 2, size=(n_samples, n_labels))
y_pred = bernoulli.rvs(np.ones((n_samples, n_labels)) / 2, size=(n_samples, n_labels))

@profile
def wrap_precision_recall_fscore_support(y_true, y_pred):
    return precision_recall_fscore_support(y_true, y_pred)

@profile
def wrap_precision_recall_fscore_support_with_multilabel_confusion_matrix(y_true, y_pred):
    return precision_recall_fscore_support_with_multilabel_confusion_matrix(y_true, y_pred)

wrap_precision_recall_fscore_support(y_true, y_pred)
wrap_precision_recall_fscore_support_with_multilabel_confusion_matrix(y_true, y_pred)
