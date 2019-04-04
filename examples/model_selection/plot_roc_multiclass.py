"""
==================================================
Multiclass Receiver Operating Characteristic (ROC)
==================================================

The multiclass One-vs-Rest scheme is functionally the same as the multilabel
case with one label per sample. See
:ref:`sphx_glr_auto_examples_model_selection_plot_roc.py` for an example
of multilabel roc.

The mutliclass One-vs-One scheme compares every unique pairwise combination
of classes. In this example, we train a SVM classifier on the iris dataset and
calcuate the AUC using a macro average and a average weighted by the number of
true instances for each label combination.

.. topic:: References:

    .. [HT2001] Hand, D.J. and Till, R.J., 2001. `A simple generalisation
       of the area under the ROC curve for multiple class classification
       problems. <http://link.springer.com/article/10.1023/A:1010920819831>`_
       Machine learning, 45(2), pp.171-186.
"""
print(__doc__)

###############################################################################
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
y_prob = classifier.fit(X_train, y_train).predict_proba(X_test)

###############################################################################
# Compute the AUC score
# .....................
# The ROC area can be calculated by taking the a macro average or an
# average weighted by the number of true instances for each label combination.
from sklearn.metrics import roc_auc_score

macro_roc_auc_ovo = roc_auc_score(y_test, y_prob, multiclass="ovo",
                                  average="macro")
weighted_roc_auc_ovo = roc_auc_score(y_test, y_prob, multiclass="ovo",
                                     average="weighted")
print("One-vs-One ROC AUC scores: {0} (uniform), {1} (weighted by prevalence)"
      .format(macro_roc_auc_ovo, weighted_roc_auc_ovo))

###############################################################################
# Plot ROC curves for the multiclass problem using One-vs-One
# ...........................................................
# The ROC curve for every pair of classes are drawn together with the
# average weighted uniformly and weighted by prevalence.
import matplotlib.pyplot as plt
from itertools import combinations, permutations
from sklearn.metrics import roc_curve, auc
from scipy import interp

n_classes = len(np.unique(y))

fpr = dict()
tpr = dict()
roc_auc = dict()
prevalence = dict()
for a, b in combinations(range(n_classes), 2):
    ab_mask = np.logical_or(y_test == a, y_test == b)

    # Compute ROC curve and ROC area with `a` as the positive class
    fpr[(a, b)], tpr[(a, b)], _ = roc_curve(
          y_test[ab_mask] == a, y_prob[ab_mask, a])
    roc_auc[(a, b)] = auc(fpr[(a, b)], tpr[(a, b)])
    prevalence[(a, b)] = np.average(ab_mask)

    # Compute ROC curve and ROC area with `b` as the positive class
    fpr[(b, a)], tpr[(b, a)], _ = roc_curve(
          y_test[ab_mask] == b, y_prob[ab_mask, b])
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
          lw=2, label='class {0} vs. {1} '
          '(area = {2:0.2f})'.format(a, b, roc_auc[(a, b)]))
ax.plot([0, 1], [0, 1], 'k--', lw=2)
ax.set_xlim([0.0, 1.0])
ax.set_ylim([0.0, 1.05])
ax.set_xlabel('False Positive Rate')
ax.set_ylabel('True Positive Rate')
ax.set_title('An extension of ROC to multi-class using One-vs-One')
ax.legend(loc="lower right")
fig.show()
