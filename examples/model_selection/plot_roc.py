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

ROC curves are typically used in binary classification to study the output of
a classifier. The ROC curve and ROC area can be extended to multi-class or
multi-label classification by using the One-vs-Rest or One-vs-One scheme.

.. note::

    See also :func:`sklearn.metrics.roc_auc_score`,
             :ref:`sphx_glr_auto_examples_model_selection_plot_roc_crossval.py`.
"""
print(__doc__)
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
# The ROC area can be calculated by taking the average either unweighted
# or weighted by the number of true instances for each label.
from sklearn.metrics import roc_auc_score

y_score_norm = np.exp(y_score)/np.exp(y_score).sum(axis=-1, keepdims=True)
unweighted_roc_auc_ovr = roc_auc_score(y_test, y_score_norm, multiclass="ovr")
weighted_roc_auc_ovr = roc_auc_score(
      y_test, y_score_norm, multiclass="ovr", average="weighted")
print("One-vs-Rest ROC AUC scores: {0} (unweighted), {1} (weighted)".format(
       unweighted_roc_auc_ovr, weighted_roc_auc_ovr))

###############################################################################
# Plotting the ROC curve for virginica
# ....................................
# A ROC curve is drawn by considering each element of the label indicator
# matrix as a binary prediction (micro-averaging). In the following, the ROC
# curve for virginica is drawn.
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

fig, ax = plt.subplots()
lw = 2
ax.plot(fpr[2], tpr[2], color='darkorange',
        lw=lw, label='ROC curve (area = %0.2f)' % roc_auc[2])
ax.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
ax.set_xlim([0.0, 1.0])
ax.set_ylim([0.0, 1.05])
ax.set_xlabel('False Positive Rate')
ax.set_ylabel('True Positive Rate')
ax.set_title('Receiver operating characteristic example for {}'.format(
    iris.target_names[2]))
ax.legend(loc="lower right")
fig.show()

###############################################################################
# Plot ROC curves for the multiclass problem using One-vs-Rest
# ............................................................
# A ROC curve is drawn using macro-averaging, which gives equal weight to the
# classification of each label.
from itertools import cycle
from scipy import interp

# First aggregate all false positive rates
all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))

# Then interpolate all ROC curves at these points
mean_tpr = np.zeros_like(all_fpr)
for i in range(n_classes):
    mean_tpr += interp(all_fpr, fpr[i], tpr[i])

# Finally average it and compute AUC
mean_tpr /= n_classes

fpr["macro"] = all_fpr
tpr["macro"] = mean_tpr
roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

# Plot all ROC curves
fig, ax = plt.subplots()
ax.plot(fpr["micro"], tpr["micro"],
        label='micro-average ROC curve (area = {0:0.2f})'
              ''.format(roc_auc["micro"]),
        color='deeppink', linestyle=':', linewidth=4)

ax.plot(fpr["macro"], tpr["macro"],
        label='macro-average ROC curve (area = {0:0.2f})'
              ''.format(roc_auc["macro"]),
        color='navy', linestyle=':', linewidth=4)

colors = cycle(['aqua', 'darkorange', 'cornflowerblue'])
for i, color in zip(range(n_classes), colors):
    plt.plot(fpr[i], tpr[i], color=color, lw=lw,
             label='ROC curve of class {0} (area = {1:0.2f})'
             ''.format(i, roc_auc[i]))

ax.plot([0, 1], [0, 1], 'k--', lw=lw)
ax.set_xlim([0.0, 1.0])
ax.set_ylim([0.0, 1.05])
ax.set_xlabel('False Positive Rate')
ax.set_ylabel('True Positive Rate')
ax.set_title('An extension of ROC to multi-class using One-vs-Rest')
ax.legend(loc="lower right")
fig.show()

###############################################################################
# One-vs-One
# ---------------------
# The One-vs-One scheme compares every unique pairwise combination of classes.
#
# Compute the AUC score
# .....................
# The ROC area can be calculated by taking the average either unweighted
# or weighted by the number of true instances for each label.
macro_roc_auc_ovo = roc_auc_score(
      y_test, y_score_norm, multiclass="ovo", average="macro")
weighted_roc_auc_ovo = roc_auc_score(
      y_test, y_score_norm, multiclass="ovo", average="weighted")
print("One-vs-One ROC AUC scores: {0} (uniform), {1} (weighted by prevalence)"
      .format(macro_roc_auc_ovo, weighted_roc_auc_ovo))

###############################################################################
# Plot ROC curves for the multiclass problem using One-vs-One
# ...........................................................
# The ROC curve for every pair of classes are drawn together with the
# average weighted uniformly and weighted by prevalence.
from itertools import combinations, permutations

prevalence = dict()
for a, b in combinations(range(n_classes), 2):
    ab_mask = np.logical_or(y_test == a, y_test == b)

    # Compute ROC curve and ROC area with `a` as the positive class
    fpr[(a, b)], tpr[(a, b)], _ = roc_curve(
          y_test[ab_mask] == a, y_score[ab_mask, a])
    roc_auc[(a, b)] = auc(fpr[(a, b)], tpr[(a, b)])
    prevalence[(a, b)] = np.average(ab_mask)

    # Compute ROC curve and ROC area with `b` as the positive class
    fpr[(b, a)], tpr[(b, a)], _ = roc_curve(
          y_test[ab_mask] == b, y_score[ab_mask, b])
    roc_auc[(b, a)] = auc(fpr[(b, a)], tpr[(b, a)])
    prevalence[(b, a)] = np.average(ab_mask)

class_permutations = list(permutations(range(n_classes), 2))
all_multiclass_fpr = np.unique(
      np.concatenate([fpr[(a, b)] for a, b in class_permutations]))

multiclass_interp_tpr = dict()
for a, b in class_permutations:
    multiclass_interp_tpr[(a, b)] = interp(
          all_multiclass_fpr, fpr[(a, b)], tpr[(a, b)])

all_multiclass_tpr = np.array(
      [multiclass_interp_tpr[(a, b)] for a, b in class_permutations])
all_prevalence = np.array([prevalence[(a, b)] for a, b in class_permutations])

roc_auc_uniform_average_tpr = np.average(all_multiclass_tpr, axis=0)
roc_auc_prevalence_average_tpr = np.average(
      all_multiclass_tpr, axis=0, weights=all_prevalence)


fig, ax = plt.subplots()
# plot roc curve as a macro average
ax.plot(
    all_multiclass_fpr,
    roc_auc_uniform_average_tpr,
    color='navy',
    linestyle=':',
    lw=4,
    label='macro average (area = {0:0.2f})'.format(
        macro_roc_auc_ovo),
)
# plot roc curve as a weighted average
ax.plot(
    all_multiclass_fpr,
    roc_auc_prevalence_average_tpr,
    color='deeppink',
    linestyle=':',
    lw=4,
    label='weighted average (area = {0:0.2f})'.format(
        weighted_roc_auc_ovo),
)

# plot roc curve for every of classes
for a, b in permutations(range(n_classes), 2):
    ax.plot(
          fpr[(a, b)], tpr[(a, b)],
          lw=lw, label='class {0} vs. {1} '
          '(area = {2:0.2f})'.format(a, b, roc_auc[(a, b)]))
ax.plot([0, 1], [0, 1], 'k--', lw=lw)
ax.set_xlim([0.0, 1.0])
ax.set_ylim([0.0, 1.05])
ax.set_xlabel('False Positive Rate')
ax.set_ylabel('True Positive Rate')
ax.set_title('An extension of ROC to multi-class using One-vs-One')
ax.legend(loc="lower right")
fig.show()
