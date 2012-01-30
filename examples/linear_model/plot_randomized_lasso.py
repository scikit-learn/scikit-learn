"""
==============================================
Randomized Lasso: feature selection with Lasso
==============================================

Performs feature scoring and selection using a :ref:`randomized sparse
linear model <randomized_l1>` (Lasso).

The problem here is, given a small number of observations, to recover
which features of X are relevant to explain y. For this we use the
randomized lasso, which can outperform standard statistical tests if the
true model is sparse, i.e. if a small fraction of the features are
relevant. The output of this procedure is a stability score for each
feature: the higher, the more likely the feature is to be relevant.

The use of the randomized lasso requires the choice of the alpha
parameter setting the sparsity of the model estimated. In the first
figure, we vary this parameter and look at the stability score as a
function of it. This analysis, knowing the ground truth shows an optimal
regime in which relevant features stand out from the irrelevant ones.

Here, we set alpha by cross-validation. Although this choice may not lead
to a optimal choice to separate relevant from irrelevant features, the
stability score of the randomized lasso still outperforms standard
statistics, as shown on the second figure. To quantify the performance of
different feature selection methods, with use the area under curve (AUC)
of the precision-recall.

Note, without knowing the ground truth, setting the value of alpha is
challenging. Choosing alpha by cross-validation to minimize left-out
residuals leads to under-penalized models: including a small number of
non-relevant variables is not detrimental to prediction score. For very
noisy data, or very few samples, it may be useful to use AIC to set
alpha, which tends, on the opposite, to set high values of alpha.
However, if alpha is too high, the LARS is similar to a step-wise
regression, and it will not display any gain over univariate
feature selection.

XXX: discuss the parameters in the choice of the simulation

* Sparsity, n_samples, n_features
* Conditionning of the design matrix on the signal subspace
* Amount of noise
"""
#print __doc__

import pylab as pl
import numpy as np
from scipy import linalg

from sklearn.linear_model import RandomizedLasso, lasso_stability_path, \
                                 LassoLarsCV
from sklearn.cross_validation import ShuffleSplit
from sklearn.feature_selection import f_regression
from sklearn.preprocessing import Scaler
from sklearn.metrics import auc, precision_recall_curve

###############################################################################
# Simulate regression data with a correlated design
n_features = 201 # Use a multiple of 3
n_relevant_features = 10
noise_level = .1
coef_min = 1
# The Donoho-Tanner phase transition is around n_samples=35: below we
# will completely fail to recover in the well-conditionned case
n_samples = 40
conditionning = 1e-4
block_size = 5

rng = np.random.RandomState(42)

# The coefficients of our model
coef = np.zeros(n_features)
coef[:n_relevant_features] = coef_min + rng.rand(n_relevant_features)

# The correlation of our design: variables correlated by blocs of 3
corr = np.zeros((n_features, n_features))
for i in range(0, n_features, block_size):
    corr[i:i + block_size, i:i + block_size] = 1 - conditionning
corr.flat[::n_features + 1] = 1
corr = linalg.cholesky(corr)

# Our design
X = rng.normal(size=(n_samples, n_features))
X = np.dot(X, corr)
X = Scaler().fit_transform(X.copy())

# The output variable
y = np.dot(X, coef)
y -= np.mean(y)
y /= np.std(y)
# We scale the added noise as a function of the average correlation
# between the design and the output variable
y += noise_level * rng.normal(size=n_samples)


###############################################################################
# Plot stability selection path
alpha_grid, scores_path = lasso_stability_path(X, y,
                                               random_state=42, eps=0.01)

pl.figure()
# We plot the path as a function of alpha/alpha_max to the power 1/3: the
# power 1/3 scales the path less brutally than the log, and enables to
# see the progression along the path
hg = pl.plot(alpha_grid[1:] ** .333, scores_path[coef != 0].T[1:], 'r')
hb = pl.plot(alpha_grid[1:] ** .333, scores_path[coef == 0].T[1:], 'k')
ymin, ymax = pl.ylim()
pl.xlabel(r'$(\alpha / \alpha_{max})^{1/3}$')
pl.ylabel('Stability score: proportion of times selected')
pl.title('Stability Scores Path')
pl.axis('tight')
pl.legend((hg[0], hb[0]), ('relevant features', 'irrelevant features'),
          loc='best')
pl.title('Conditionning %e' % conditionning)

###############################################################################
# Plot the estimated stability scores for best cross-validated alpha

# First find the best alpha:
cv = ShuffleSplit(n_samples, test_fraction=.25, random_state=42)
lars_cv = LassoLarsCV(cv=cv).fit(X, y)
alpha_2 = lars_cv.alphas_[0]
alpha_1 = .1*alpha_2
alphas = np.linspace(alpha_1, alpha_2, 6)
# Then run the RandomizedLasso
clf = RandomizedLasso(alpha=alphas, random_state=42)
clf.fit(X, y)

F, _ = f_regression(X, y)  # compare with F-score

pl.figure()
for name, score in [('F-score', F), ('Stability selection (max)', clf.scores_),
                ('Stability selection (mean)', clf.all_scores_.mean(axis=-1)),
                ('L1 coefs', np.abs(lars_cv.coef_))]:
    precision, recall, thresholds = precision_recall_curve(coef != 0,
                score)
    area = auc(recall, precision)
    pl.semilogy(np.maximum(score / np.max(score), 1e-4),
            label="%s. AUC: %.3f" % (name, area))

pl.plot(np.where(coef != 0)[0], [2e-4] * n_relevant_features, 'mo',
        label="Ground truth")
pl.xlabel("Features")
pl.ylabel("Score")
# Plot only the 100 first coefficients
pl.xlim(0, 100)
pl.legend(loc='lower right')
pl.title('Conditionning %e' % conditionning)
pl.show()

