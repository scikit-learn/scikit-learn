'''
==================
RBF SVM parameters
==================

This example illustrates the effect of the parameters `gamma`
and `C` of the rbf kernel SVM.

Intuitively, the `gamma` parameter defines how far the influence
of a single training example reaches, with low values meaning 'far'
and high values meaning 'close'.
The `C` parameter trades off misclassification of training examples
against simplicity of the decision surface. A low C makes
the decision surface smooth, while a high C aims at classifying
all training examples correctly.

Two plots are generated. The first is a visualization of the decision function
for a variety of parameter values, and the second is a heatmap of the
classifier's cross-validation accuracy and training time as a function of `C`
and `gamma`.

An interesting observation on overfitting can be made when comparing validation
and training error: higher C always result in lower training error, as it
inceases complexity of the classifier.

For the validation set on the other hand, there is a tradeoff of goodness of
fit and generalization.

We can observe that the lower right half of the parameters (below the diagonal
with high C and gamma values) is characteristic of parameters that yields an
overfitting model: the trainin score is very high but there is a wide gap. The
top and left parts of the parameter plots show underfitting models: the C and
gamma values can individually or in conjunction constrain the model too much
leading to low training scores (hence low validation scores too as validation
scores are on average upper bounded by training scores).


We can also see that the training time is quite sensitive to the parameter
setting, while the prediction time is not impacted very much. This is probably
a consequence of the small size of the data set.
'''
print(__doc__)

import numpy as np
import pylab as pl

from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import load_iris
from sklearn.cross_validation import StratifiedKFold
from sklearn.grid_search import GridSearchCV

##############################################################################
# Load and prepare data set
#
# dataset for grid search
iris = load_iris()
X = iris.data
Y = iris.target

# dataset for decision function visualization
X_2d = X[:, :2]
X_2d = X_2d[Y > 0]
Y_2d = Y[Y > 0]
Y_2d -= 1

# It is usually a good idea to scale the data for SVM training.
# We are cheating a bit in this example in scaling all of the data,
# instead of fitting the transformation on the training set and
# just applying it on the test set.

scaler = StandardScaler()

X = scaler.fit_transform(X)
X_2d = scaler.fit_transform(X_2d)

##############################################################################
# Train classifier
#
# For an initial search, a logarithmic grid with basis
# 10 is often helpful. Using a basis of 2, a finer
# tuning can be achieved but at a much higher cost.

C_range = 10.0 ** np.arange(-2, 9)
gamma_range = 10.0 ** np.arange(-5, 4)
param_grid = dict(gamma=gamma_range, C=C_range)
cv = StratifiedKFold(y=Y, n_folds=3)
grid = GridSearchCV(SVC(), param_grid=param_grid, cv=cv,
                    compute_training_score=True)
grid.fit(X, Y)

print("The best classifier is: ", grid.best_estimator_)

# Now we need to fit a classifier for all parameters in the 2d version
# (we use a smaller set of parameters here because it takes a while to train)
C_2d_range = [1, 1e2, 1e4]
gamma_2d_range = [1e-1, 1, 1e1]
classifiers = []
for C in C_2d_range:
    for gamma in gamma_2d_range:
        clf = SVC(C=C, gamma=gamma)
        clf.fit(X_2d, Y_2d)
        classifiers.append((C, gamma, clf))

##############################################################################
# visualization
#
# draw visualization of parameter effects
pl.figure(figsize=(8, 6))
xx, yy = np.meshgrid(np.linspace(-5, 5, 200), np.linspace(-5, 5, 200))
for (k, (C, gamma, clf)) in enumerate(classifiers):
    # evaluate decision function in a grid
    Z = clf.decision_function(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)

    # visualize decision function for these parameters
    pl.subplot(len(C_2d_range), len(gamma_2d_range), k + 1)
    pl.title("gamma 10^%d, C 10^%d" % (np.log10(gamma), np.log10(C)),
             size='medium')

    # visualize parameter's effect on decision function
    pl.pcolormesh(xx, yy, -Z, cmap=pl.cm.jet)
    pl.scatter(X_2d[:, 0], X_2d[:, 1], c=Y_2d, cmap=pl.cm.jet)
    pl.xticks(())
    pl.yticks(())
    pl.axis('tight')

# plot the scores of the grid
# cv_scores_ contains parameter settings and scores
score_dict = grid.cv_scores_

# We extract validation and training scores, as well as training and prediction
# times
_, val_scores, _, train_scores, train_time, pred_time = zip(*score_dict)

arrays = [val_scores, train_scores, train_time, pred_time]
titles = ["Validation Score", "Training Score", "Training Time",
          "Prediction Time"]

# for each value draw heatmap  as a function of gamma and C
pl.figure(figsize=(12, 8))
for i, (arr, title) in enumerate(zip(arrays, titles)):
    pl.subplot(2, 2, i + 1)
    arr = np.array(arr).reshape(len(C_range), len(gamma_range))
    pl.title(title)
    pl.imshow(arr, interpolation='nearest', cmap=pl.cm.spectral)
    pl.xlabel('gamma')
    pl.ylabel('C')
    pl.colorbar()
    pl.xticks(np.arange(len(gamma_range)), ["%.e" % g for g in gamma_range],
              rotation=45)
    pl.yticks(np.arange(len(C_range)), ["%.e" % C for C in C_range])

pl.subplots_adjust(top=.95, hspace=.35, left=.0, right=.8, wspace=.05)

pl.show()
