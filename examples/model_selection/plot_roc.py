"""
=======================================
Receiver Operating Characteristic (ROC)
=======================================

Example of Receiver Operating Characteristic (ROC) metric to evaluate
classifier output quality.

ROC curves typically feature true positive rate on the Y axis, and false
positive rate on the X axis. This means that the top left corner of the plot is
the "ideal" point - a false positive rate of zero, and a true positive rate of
one. This is not very realistic, but it does mean that a larger area under the
curve (AUC) is usually better.

The "steepness" of ROC curves is also important, since it is ideal to maximize
the true positive rate while minimizing the false positive rate.

Multiclass settings
-------------------

ROC curves are typically used in binary classification to study the output of
a classifier. Extensions of ROC curve and ROC area to multi-class
or multi-label classification can use the One-vs-Rest or One-vs-One scheme.

One-vs-Rest
-----------

The output is binarized and one ROC curve can be drawn per label,
where the label is the positive class and all other labels are
the negative class.

The ROC area can be approximated by taking the average--unweighted or weighted
by the a priori class distribution--of the one-vs-rest ROC areas.

One can also draw a ROC curve by considering
each element of the label indicator matrix as a binary prediction
(micro-averaging).

Another evaluation measure for one-vs-rest multi-class classification is
macro-averaging, which gives equal weight to the classification of each
label.

One-vs-One
----------

Two ROC curves can be drawn per pair of labels because either of the two
labels can be considered the positive class.

The ROC area can be approximated by first computing the
approximate ROC area of each label pair as the average of the
two ROC AUC scores corresponding to that pair. The One-vs-One
approximation of a multi-class ROC AUC score is the average--
unweighted or weighted by the a priori class distribution--across
all of the pairwise approximate ROC AUC scores.

.. note::

    See also :func:`sklearn.metrics.roc_auc_score`,
             :ref:`sphx_glr_auto_examples_model_selection_plot_roc_crossval.py`.

"""
print(__doc__)

import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle

from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc, roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from scipy import interp

# Import some data to play with
iris = datasets.load_iris()
X = iris.data
y = iris.target

classes = np.unique(y)
n_classes = len(classes)

# Add noisy features to make the problem harder
random_state = np.random.RandomState(0)
n_samples, n_features = X.shape
X = np.c_[X, random_state.randn(n_samples, 200 * n_features)]

# shuffle and split training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.5,
                                                    random_state=0)

# Learn to predict each class against the other
classifier = OneVsRestClassifier(svm.SVC(kernel='linear', probability=True,
                                 random_state=random_state))
y_score = classifier.fit(X_train, y_train).decision_function(X_test)

# Compute ROC curve and ROC area for each class

# Binarize y_test to compute the ROC curve
y_test_binarized = label_binarize(y_test, classes=classes)

fpr = dict()
tpr = dict()
roc_auc = dict()
for i in range(n_classes):
    fpr[i], tpr[i], _ = roc_curve(y_test_binarized[:, i], y_score[:, i])
    roc_auc[i] = auc(fpr[i], tpr[i])


##############################################################################
# Plot of a ROC curve for a specific class
plt.figure()
lw = 2
plt.plot(fpr[2], tpr[2], color='darkorange',
         lw=lw, label='ROC curve (area = %0.2f)' % roc_auc[2])
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example')
plt.legend(loc="lower right")
plt.show()


##############################################################################
# Plot ROC curves for the multiclass problem using One vs. Rest classification.

# Compute micro-average ROC curve and ROC area
fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel())
roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

# Compute macro-average ROC curve and ROC area

# First aggregate all false positive rates
all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))

# Then interpolate all ROC curves at this points
mean_tpr = np.zeros_like(all_fpr)
for i in range(n_classes):
    mean_tpr += interp(all_fpr, fpr[i], tpr[i])

# Finally average it and compute AUC
mean_tpr /= n_classes

fpr["macro"] = all_fpr
tpr["macro"] = mean_tpr
roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

# Plot all ROC curves
plt.figure()
plt.plot(fpr["micro"], tpr["micro"],
         label='micro-average ROC curve (area = {0:0.2f})'
               ''.format(roc_auc["micro"]),
         color='deeppink', linestyle=':', linewidth=4)

plt.plot(fpr["macro"], tpr["macro"],
         label='macro-average ROC curve (area = {0:0.2f})'
               ''.format(roc_auc["macro"]),
         color='navy', linestyle=':', linewidth=4)

colors = cycle(['aqua', 'darkorange', 'cornflowerblue'])
for i, color in zip(range(n_classes), colors):
    plt.plot(fpr[i], tpr[i], color=color, lw=lw,
             label='ROC curve of class {0} (area = {1:0.2f})'
             ''.format(i, roc_auc[i]))

plt.plot([0, 1], [0, 1], 'k--', lw=lw)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('An extension of Receiver operating characteristic to multi-class '
          'using One-vs-Rest')
plt.legend(loc="lower right")
plt.show()

# Compute the One-vs-Rest ROC AUC score, weighted and unweighted
unweighted_roc_auc_ovr = roc_auc_score(y_test, y_score, multiclass="ovr")
weighted_roc_auc_ovr = roc_auc_score(
    y_test, y_score, multiclass="ovr", average="weighted")
print("One-vs-Rest ROC AUC scores: {0} (unweighted), {1} (weighted)".format(
      unweighted_roc_auc_ovr, weighted_roc_auc_ovr))

##############################################################################
# Plot ROC curves for the multiclass problem using One vs. One classification.

for pos in range(n_classes):
    for neg in range(pos + 1, n_classes):
        # Filter `y_test` and `y_score` to only consider the current
        # class pair: `pos` and `neg`.
        class_pair_indices = np.in1d(y_test, [pos, neg])
        y_true_filtered = y_test[class_pair_indices]
        y_score_filtered = y_score[class_pair_indices]

        # Compute ROC curve and ROC area with `pos` as the positive class
        class_a = y_true_filtered == pos
        fpr[(pos, neg)], tpr[(pos, neg)], _ = roc_curve(
                class_a, y_score_filtered[:, pos])
        roc_auc[(pos, neg)] = auc(fpr[(pos, neg)], tpr[(pos, neg)])

        # Compute ROC curve and ROC area with `neg` as the positive class
        class_b = y_true_filtered == neg
        fpr[(neg, pos)], tpr[(neg, pos)], _ = roc_curve(
                class_b, y_score_filtered[:, neg])
        roc_auc[(neg, pos)] = auc(fpr[(neg, pos)], tpr[(neg, pos)])

plt.figure()
for pos in range(n_classes):
    for neg in range(pos + 1, n_classes):
        plt.plot(fpr[(pos, neg)], tpr[(pos, neg)], lw=lw,
                 label='ROC curve of class {0} against class {1} '
                       '(area = {2:0.2f})'.format(
                        pos, neg, roc_auc[(pos, neg)]))
        plt.plot(fpr[(neg, pos)], tpr[(neg, pos)], lw=lw,
                 label='ROC curve of class {0} against class {1} '
                       '(area = {2:0.2f})'.format(
                        neg, pos, roc_auc[(neg, pos)]))
plt.plot([0, 1], [0, 1], 'k--', lw=lw)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('An extension of Receiver operating characteristic to multi-class '
          'using One-vs-One')
plt.legend(bbox_to_anchor=(1.8, 0.55))
plt.show()

# Compute the One-vs-One ROC AUC score, weighted and unweighted
unweighted_roc_auc_ovo = roc_auc_score(y_test, y_score, multiclass="ovo")
weighted_roc_auc_ovo = roc_auc_score(
    y_test, y_score, multiclass="ovo", average="weighted")
print("One-vs-One ROC AUC scores: {0} (unweighted), {1} (weighted)".format(
      unweighted_roc_auc_ovo, weighted_roc_auc_ovo))
