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
a classifier. The ROC curve and ROC area can be extended to multi-class or
multi-label classification by using the One-vs-Rest or One-vs-One scheme.

.. note::

    See also :func:`sklearn.metrics.roc_auc_score`,
             :ref:`sphx_glr_auto_examples_model_selection_plot_roc_crossval.py`.
"""
###############################################################################
# One-vs-Rest
# -----------
# The One-vs-Rest scheme compares the each class against all of the other
# classes ("the rest").
#
# Load iris dataset and train a SVC
# .................................
from sklearn import svm, datasets
from sklearn.model_selection import train_test_split
from sklearn.multiclass import OneVsRestClassifier
import numpy as np

# Import some data to play with
iris = datasets.load_iris()
X = iris.data
y = iris.target

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

###############################################################################
# Compute the AUC scores
# ......................
# The ROC area can be approximated by taking the average either weighted
# uniformly or by the priori class distribution.
from sklearn.metrics import roc_auc_score

y_score_norm = y_score / y_score.sum(1, keepdims=True)
unweighted_roc_auc_ovr = roc_auc_score(y_test, y_score_norm, multiclass="ovr")
weighted_roc_auc_ovr = roc_auc_score(
      y_test, y_score_norm, multiclass="ovr", average="weighted")
print("One-vs-Rest ROC AUC scores: {0} (unweighted), {1} (weighted)".format(
       unweighted_roc_auc_ovr, weighted_roc_auc_ovr))

###############################################################################
# Plotting the ROC curve for virginica
# ....................................
import matplotlib.pyplot as plt
from sklearn.preprocessing import label_binarize
from sklearn.metrics import roc_curve, auc

# Binarize y_test to compute the ROC curve
classes = np.unique(y)
n_classes = len(classes)
y_test_binarized = label_binarize(y_test, classes=classes)

# Compute ROC curve and ROC area for each class
fpr = dict()
tpr = dict()
roc_auc = dict()
for i in range(n_classes):
    fpr[i], tpr[i], _ = roc_curve(y_test_binarized[:, i], y_score[:, i])
    roc_auc[i] = auc(fpr[i], tpr[i])

fpr["micro"], tpr["micro"], _ = roc_curve(
      y_test_binarized.ravel(), y_score.ravel())
roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

plt.figure()
lw = 2
plt.plot(fpr[2], tpr[2], color='darkorange',
         lw=lw, label='ROC curve (area = %0.2f)' % roc_auc[2])
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example for {}'
          .format(iris.target_names[2]))
plt.legend(loc="lower right")

###############################################################################
# Plot ROC curves for the multiclass problem using One-vs-Rest
# ............................................................
# One can draw a ROC curve by considering each element of the label indicator
# matrix as a binary prediction (micro-averaging).
#
# Another evaluation measure for one-vs-rest multi-class classification is
# macro-averaging, which gives equal weight to the classification of each
# label.
from itertools import cycle
from scipy import interp

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
plt.title('An extension of ROC to multi-class using One-vs-Rest')
plt.legend(loc="lower right")

###############################################################################
# One-vs-One
# ---------------------
# The One-vs-One scheme compares every pairwise combiantion of classes.
#
# Compute the AUC score
# .....................
# The ROC area can be approximated by taking the average either weighted
# uniformly or by the priori class distribution.
unweighted_roc_auc_ovo = roc_auc_score(y_test, y_score_norm, multiclass="ovo")
weighted_roc_auc_ovo = roc_auc_score(
      y_test, y_score_norm, multiclass="ovo", average="weighted")
print("One-vs-One ROC AUC scores: {0} (unweighted), {1} (weighted)".format(
       unweighted_roc_auc_ovo, weighted_roc_auc_ovo))

###############################################################################
# Plot ROC curves for the multiclass problem using One-vs-One
# ...........................................................
from itertools import combinations

for a, b in combinations(range(n_classes), 2):
    ab_mask = np.logical_or(y_test == a, y_test == b)

    # Compute ROC curve and ROC area with `a` as the positive class
    fpr[(a, b)], tpr[(a, b)], _ = roc_curve(
          y_test[ab_mask] == a, y_score[ab_mask, a])
    roc_auc[(a, b)] = auc(fpr[(a, b)], tpr[(a, b)])

    # Compute ROC curve and ROC area with `b` as the positive class
    fpr[(b, a)], tpr[(b, a)], _ = roc_curve(
          y_test[ab_mask] == b, y_score[ab_mask, b])
    roc_auc[(b, a)] = auc(fpr[(b, a)], tpr[(b, a)])

plt.figure()
for a, b in combinations(range(n_classes), 2):
    plt.plot(
          fpr[(a, b)], tpr[(a, b)],
          lw=lw, label='ROC curve: class {0} vs. {1} '
          '(area = {2:0.2f})'.format(a, b, roc_auc[(a, b)]))
    plt.plot(
           fpr[(b, a)], tpr[(b, a)],
           lw=lw, label='ROC curve: class {0} vs. {1} '
           '(area = {2:0.2f})'.format(b, a, roc_auc[(b, a)]))
plt.plot([0, 1], [0, 1], 'k--', lw=lw)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('An extension of ROC to multi-class using One-vs-One')
plt.legend(loc="lower right")
