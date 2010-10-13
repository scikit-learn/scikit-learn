# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Mathieu Blondel <mathieu@mblondel.org>
# License: BSD Style.

import numpy as np

def confusion_matrix(y, y_):
    """
    compute confusion matrix
    to evaluate the accuracy of a classification result

    By definition a confusion matrix cm is such that

    cm[i,j] is equal to the number of observations known to be in group i
    but predicted to be in group j

    Parameters
    ==========

    y : array, shape = [n_samples]
        true targets

    y_ : array, shape = [n_samples]
        estimated targets

    Returns
    =======
    cm : array, shape = [n_classes,n_classes]
        confusion matrix

    References
    ==========
    http://en.wikipedia.org/wiki/Confusion_matrix
    """
    # removing possible NaNs in targets (they are ignored)
    clean_y = y[np.isfinite(y)].ravel()
    clean_y_ = y_[np.isfinite(y_)].ravel()

    labels = np.r_[np.unique(clean_y).ravel(),np.unique(clean_y_).ravel()]
    labels = np.unique(labels)
    n_labels = labels.size

    cm = np.empty((n_labels,n_labels))
    for i, label_i in enumerate(labels):
        for j, label_j in enumerate(labels):
            cm[i,j] = np.sum(np.logical_and(y==label_i, y_==label_j))

    return cm


def roc(y, probas_):
    """compute Receiver operating characteristic (ROC)

    Parameters
    ==========

    y : array, shape = [n_samples]
        true targets

    probas_ : array, shape = [n_samples]
        estimated probabilities

    Returns
    =======
    fpr : array, shape = [n]
        False Positive Rates

    tpr : array, shape = [n]
        True Positive Rates

    thresholds : array, shape = [n]
        Thresholds on proba_ used to compute fpr and tpr

    References
    ==========
    http://en.wikipedia.org/wiki/Receiver_operating_characteristic
    """
    y = y.ravel()
    probas_ = probas_.ravel()
    thresholds = np.sort(np.unique(probas_))[::-1]
    n_thresholds = thresholds.size
    tpr = np.empty(n_thresholds) # True positive rate
    fpr = np.empty(n_thresholds) # False positive rate
    n_pos = float(np.sum(y==1)) # nb of true positive
    n_neg = float(np.sum(y==0)) # nb of true negative
    for i, t in enumerate(thresholds):
        tpr[i] = np.sum(y[probas_>=t]==1) / n_pos
        fpr[i] = np.sum(y[probas_>=t]==0) / n_neg

    return fpr, tpr, thresholds


def auc(x, y):
    """Compute Area Under the Curve (AUC)
    using the trapezoidal rule

    Parameters
    ==========

    x : array, shape = [n]
        x coordinates

    y : array, shape = [n]
        y coordinates

    Returns
    =======
    auc : float

    """
    x = np.asanyarray(x)
    y = np.asanyarray(y)
    h = np.diff(x)
    area = np.sum(h * (y[1:]+y[:-1])) / 2.0
    return area

def precision(y_true, y_pred):
    """
    Compute precision.

    Parameters
    ==========
    y_true : array, shape = [n_samples]
        true targets

    y_pred : array, shape = [n_samples]
        predicted targets

    Returns
    =======
    precision : float
    """
    true_pos = np.sum(y_true[y_pred == 1]==1)
    false_pos = np.sum(y_true[y_pred == 1]==0)
    return true_pos / float(true_pos + false_pos)

def recall(y_true, y_pred):
    """
    Compute recall.

    Parameters
    ==========
    y_true : array, shape = [n_samples]
        true targets

    y_pred : array, shape = [n_samples]
        predicted targets

    Returns
    =======
    recall : float
    """
    true_pos = np.sum(y_true[y_pred == 1]==1)
    false_neg = np.sum(y_true[y_pred == 0]==1)
    return true_pos / float(true_pos + false_neg)

def precision_recall(y_true, y_pred):
    """
    Compute precision and recall.

    Parameters
    ==========
    y_true : array, shape = [n_samples]
        true targets

    y_pred : array, shape = [n_samples]
        predicted targets

    Returns
    =======
    precision: float
    recall : float

    References
    ==========
    http://en.wikipedia.org/wiki/Precision_and_recall
    """
    true_pos = np.sum(y_true[y_pred == 1]==1)
    false_pos = np.sum(y_true[y_pred == 1]==0)
    false_neg = np.sum(y_true[y_pred == 0]==1)
    precision = true_pos / float(true_pos + false_pos)
    recall = true_pos / float(true_pos + false_neg)
    return precision, recall

def precision_recall_curve(y, probas_):
    """
    Compute precision-recall pairs for different probability thresholds.

    Parameters
    ==========
    y : array, shape = [n_samples]
        true targets

    probas_ : array, shape = [n_samples]
        estimated probabilities

    Returns
    =======
    precision : array, shape = [n]
        Precision values

    recall : array, shape = [n]
        Recall values

    thresholds : array, shape = [n]
        Thresholds on proba_ used to compute precision and recall
    """
    y = y.ravel()
    probas_ = probas_.ravel()
    thresholds = np.sort(np.unique(probas_))
    n_thresholds = thresholds.size + 1
    precision = np.empty(n_thresholds)
    recall = np.empty(n_thresholds)
    for i, t in enumerate(thresholds):
        y_pred = np.ones(len(y))
        y_pred[probas_ < t] = 0
        precision[i], recall[i] = precision_recall(y, y_pred)
    precision[-1] = 1.0
    recall[-1] = 0.0
    return precision, recall, thresholds

def fbeta_score(y_true, y_pred, beta):
    """
    Compute fbeta score.

    Parameters
    ==========
    y_true : array, shape = [n_samples]
        true targets

    y_pred : array, shape = [n_samples]
        predicted targets

    beta: float

    Returns
    =======
    fbeta_score: float
    """
    assert(beta > 0)
    p, r = precision_recall(y_true, y_pred)
    beta2 = beta ** 2
    return (1+beta2) * (p * r) / (beta2 * p + r)

def f1_score(y_true, y_pred):
    """
    Compute f1 score.

    Parameters
    ==========
    y_true : array, shape = [n_samples]
        true targets

    y_pred : array, shape = [n_samples]
        predicted targets

    Returns
    =======
    f1_score: float

    References
    ==========
    http://en.wikipedia.org/wiki/F1_score
    """
    return fbeta_score(y_true, y_pred, 1)


###############################################################################
# Loss functions


def zero_one(y_true, y_pred):
    """Zero-One loss
    returns the number of differences
    """
    return np.sum(y_pred != y_true)


def mean_square_error(y_true, y_pred):
    """Mean Square Error
    returns the mean square error
    """
    return np.linalg.norm(y_pred != y_true) ** 2


def explained_variance(y_true, y_pred):
    """Explained variance
    returns the explained variance
    """
    return 1 - np.var(y_true - y_pred) / np.var(y_true)

