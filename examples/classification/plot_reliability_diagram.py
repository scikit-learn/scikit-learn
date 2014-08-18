"""
=========================
Plot reliability diagram
=========================

This script generates reliability diagrams for some classifiers on an
artificial data set. Reliability diagrams allow checking if the predicted
probabilities of a binary classifier are well calibrated. For perfectly
calibrated predictions, the curve in a reliability diagram should be as close
as possible to the diagonal/identity. This would correspond to a situation in
which among $N$ instances for which a classifier predicts probability $p$ for
class $A$, the ratio of instances which actually belong to class $A$ is approx.
 $p$ (for any $p$ and sufficiently large $N$).

This script reproduces some of the results from the paper "Predicting Good
Probabilities with Supervised Learning "
http://machinelearning.wustl.edu/mlpapers/paper_files/
icml2005_Niculescu-MizilC05.pdf.

The following observations can be made:
 * Logistic regression returns well-calibrated probabilities close to the
   "perfect" line
 * Naive Bayes tends to push probabilties to 0 or 1 (note the counts in the
   histograms). This is mainly because it makes the assumption that features
   are conditionally independent given the class, which is not the case in
   this dataset which contains 2 redundant features.
 * Random Forest shows the opposite behavior: The histograms show peaks at
   approx. 0.1 and 0.9 probability, while probabilities close to 0 or 1 are
   very rare. An explanation for this is given by Niculescu-Mizil and Caruana:
   "Methods such as bagging and random forests that average predictions from
    a base set of models can have difficulty making predictions near 0
    and 1 because variance in the underlying base models will
    bias predictions that should be near zero or one away from
    these values. Because predictions are restricted to the interval [0,1],
    errors caused by variance tend to be one-sided
    near zero and one. For example, if a model should predict
    p = 0 for a case, the only way bagging can achieve this is
    if all bagged trees predict zero. If we add noise to the trees
    that bagging is averaging over, this noise will cause some
    trees to predict values larger than 0 for this case, thus moving the
    average prediction of the bagged ensemble away
    from 0. We observe this effect most strongly with random
    forests because the base-level trees trained with random
    forests have relatively high variance due to feature subseting."
   As a result, the calibration curve shows a characteristic sigmoid shape,
   indicating that the classifier could trust its "intuition" more and return
   probabilties closer to 0 or 1 typically. A post-processing such as
   Platt-calibration, which fits a sigmoid to the probabilities on a separate
   calibration dataset, would typically help if the calibration curve is
   sigmoid.
 * The scores of a Support Vector Classification (SVC), which are linearly
   related to the distance of the sample from the hyperplane, show a similar
   but even stronger effect as the Random Forest. This is not too surprising
   as the scores are in no sense probabilties and must not be interpreted as
   such as the curve shows.
 * One alternative to Platt-calibration is Isotonic Regression. While
   Platt-calibration fits a sigmoid, Isotonic Regression fits an arbitrary
   increasing (isotonic) function. Thus, it has a weaker inductive bias and
   can be applied more broadly (also in situations where the calibration curve
   is not sigmoid). The downside is that it typically requires more calibration
   data because its inductive bias is weaker. This can also be seen in the
   SVC + IR curve: While the sigmoid shape of the pure SVC scores is removed
   and the calibration curve does not show a clear bias, it is quite noisy,
   indicating much variance. Thus, the used calibration dataset (even though
   4 times larger than the training data) is too small in this case.
"""
print(__doc__)

# Author: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
# License: BSD 3 clause

import numpy as np
np.random.seed(0)
import matplotlib.pyplot as plt

from sklearn import datasets
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.isotonic import IsotonicRegression

from sklearn.metrics import reliability_curve

#### Training data
# Generate a *toy dataset* on which different classifiers are compared.
# Among the 20 features, only 2 are actually informative. 2 further features
# are redundant, i.e., linear combinations of the 2 informative features.
# The remaining features are just noise.

X, y = datasets.make_classification(n_samples=100000, n_features=20,
                                    n_informative=2, n_redundant=2)

bins = 25
# Samples used for training the models
train_samples = 100
# Additional samples for claibration using Isotonic Regression
calibration_samples = 400

X_train = X[:train_samples]
X_calibration = X[train_samples:train_samples + calibration_samples]
X_test = X[train_samples + calibration_samples:]
y_train = y[:train_samples]
y_calibration = y[train_samples:train_samples + calibration_samples]
y_test = y[train_samples + calibration_samples:]

#### Compute reliability curves for different classifiers

classifiers = {"Logistic regression": LogisticRegression(),
               "Naive Bayes": GaussianNB(),
               "Random Forest": RandomForestClassifier(n_estimators=100),
               "SVC": SVC(kernel='linear', C=1.0),
               "SVC + IR": SVC(kernel='linear', C=1.0)}

reliability_scores = {}
y_score = {}
for method, clf in classifiers.items():
    clf.fit(X_train, y_train)
    if method == "SVC + IR":  # Calibrate SVC scores using isotonic regression.
        n_plus = (y_calibration == 1.0).sum()  # number of positive examples
        n_minus = (y_calibration == 0.0).sum()  # number of negative examples
        # Determine target values for isotonic calibration. See
        # "Predicting Good Probabilities with Supervised Learning"
        # for details.
        y_target = np.where(y_calibration == 0.0,
                            1.0 / (n_minus + 2.0),
                            (n_plus + 1.0) / (n_plus + 2.0))
        # Perform actual calibration using isotonic calibration
        svm_score = clf.decision_function(X_calibration)
        ir = IsotonicRegression(out_of_bounds='clip').fit(svm_score, y_target)
        y_score[method] = ir.transform(clf.decision_function(X_test))
        reliability_scores[method] = \
            reliability_curve(y_test, y_score[method], bins=bins,
                              normalize=False)
    elif method == "SVC":
        # Use SVC scores (predict_proba returns already calibrated
        # probabilities)
        y_score[method] = clf.decision_function(X_test)
        reliability_scores[method] = \
            reliability_curve(y_test, y_score[method], bins=bins,
                              normalize=True)
    else:
        y_score[method] = clf.predict_proba(X_test)[:, 1]
        reliability_scores[method] = \
            reliability_curve(y_test, y_score[method], bins=bins,
                              normalize=False)

#### Plot reliability diagram

plt.figure(0, figsize=(8, 8))
plt.subplot2grid((3, 1), (0, 0), rowspan=2)
plt.plot([0.0, 1.0], [0.0, 1.0], 'k', label="Perfect")
for method, (bin_mean, empirical_prob_pos) in reliability_scores.items():
    scores_not_nan = np.logical_not(np.isnan(empirical_prob_pos))
    plt.plot(bin_mean[scores_not_nan],
             empirical_prob_pos[scores_not_nan], label=method)
plt.ylabel("Empirical probability")
plt.legend(loc=0)

plt.subplot2grid((3, 1), (2, 0))
for method, y_score_ in y_score.items():
    y_score_ = (y_score_ - y_score_.min()) / (y_score_.max() - y_score_.min())
    plt.hist(y_score_, range=(0, 1), bins=bins, label=method,
             histtype="step", lw=2)
plt.xlabel("Predicted Probability")
plt.ylabel("Count")
plt.legend(loc='upper center', ncol=2)
plt.show()
