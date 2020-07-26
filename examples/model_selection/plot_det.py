"""
=======================================
Detection error tradeoff (DET) curve
=======================================

In this example we compare receiver operating characteristic (ROC) and
detection error tradeoff (DET) curves to demonstrate how DET curves can help
to asses the performance of different classification algorithms.

DET curves are commonly plotted in normal deviate scale.
To achieve this we transform the errors rates as returned by the
``detection_error_tradeoff_curve`` function and the axis scale using
``scipy.stats.norm``.

.. note::

    - See :func:`sklearn.metrics.roc_curve` for further information about ROC
      curves.

    - This example is loosely based on
      :ref:`sphx_glr_auto_examples_classification_plot_classifier_comparison.py`
      .

"""
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import make_classification
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import detection_error_tradeoff_curve
from sklearn.metrics import roc_curve

from scipy.stats import norm
from matplotlib.ticker import FuncFormatter

N_SAMPLES = 1000

names = [
    "Linear SVM",
    "Random Forest",
]

classifiers = [
    SVC(kernel="linear", C=0.025),
    RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1),
]

X, y = make_classification(
    n_samples=N_SAMPLES, n_features=2, n_redundant=0, n_informative=2,
    random_state=1, n_clusters_per_class=1)

figure = plt.figure(figsize=(10, 5))

# preprocess dataset, split into training and test part
X = StandardScaler().fit_transform(X)

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=.4, random_state=0)

# prepare plots

# first prepare the ROC curve
ax_roc = plt.subplot(1, 2, 1, label='roc_curve')
ax_roc.set_title('Receiver Operating Characteristic (ROC) curves')
ax_roc.set_xlabel('False Positive Rate')
ax_roc.set_ylabel('True Positive Rate')
ax_roc.set_xlim(0, 1)
ax_roc.set_ylim(0, 1)
ax_roc.grid(linestyle='--')
ax_roc.yaxis.set_major_formatter(
    FuncFormatter(lambda y, _: '{:.0%}'.format(y)))
ax_roc.xaxis.set_major_formatter(
    FuncFormatter(lambda y, _: '{:.0%}'.format(y)))

# second prepare the DET curve
ax_det = plt.subplot(1, 2, 2, label='det_curve')
ax_det.set_title('Detection Error Tradeoff (DET) curves')
ax_det.set_xlabel('False Positive Rate')
ax_det.set_ylabel('False Negative Rate')
ax_det.set_xlim(-3, 3)
ax_det.set_ylim(-3, 3)
ax_det.grid(linestyle='--')
# customized ticks to represent normal deviate scale
ticks = [0.001, 0.01, 0.05, 0.20, 0.5, 0.80, 0.95, 0.99, 0.999]
tick_locs = norm.ppf(ticks)
tick_lbls = [
    '{:.0%}'.format(s) if (100*s).is_integer() else '{:.1%}'.format(s)
    for s in ticks
]

plt.xticks(tick_locs, tick_lbls)
plt.yticks(tick_locs, tick_lbls)

# iterate over classifiers
for name, clf in zip(names, classifiers):
    clf.fit(X_train, y_train)

    if hasattr(clf, "decision_function"):
        y_score = clf.decision_function(X_test)
    else:
        y_score = clf.predict_proba(X_test)[:, 1]

    roc_fpr, roc_tpr, _ = roc_curve(y_test, y_score)
    det_fpr, det_fnr, _ = detection_error_tradeoff_curve(y_test, y_score)

    ax_roc.plot(roc_fpr, roc_tpr)

    # transform errors into normal deviate scale
    ax_det.plot(
        norm.ppf(det_fpr),
        norm.ppf(det_fnr)
    )

# finally add legend
plt.legend(names, loc="upper right")

plt.tight_layout()
plt.show()
