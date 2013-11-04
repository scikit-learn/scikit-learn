"""
===========================================
Minimum redundancy maximum relevance (mRMR)
===========================================

Mutual information is a metric assessing the degree of statistical independence
between two random variables.

mRMR feature selection consists in selecting a subset of the available features
showing high mutual information with the target and low mutual information with
each other.

This example compares mRMR feature selection with Recursive feature elimination
(RFE) and Univariate feature selection (Uni), taking advantage of a synthetic
dataset.

This dataset has 100 samples and 3 features: A, B and C, enabling to
respectively classify 60%, 50% and 40% of the data.

Let's assume the plan is to choose only 2 of those 3 features. Given that A
and B have higher accuracy, we would expect a selection algorithm to pick those
two. However, it turns out that A and B are redundant with each other (i.e.
they are able to classify the same samples). Conversely, C has lower accuracy,
but provides indepedent information respect to A and B.

As expectable, mRMR selects feature A and C, while the other two selection
algorithm select features A and B.


.. note::

    See also :ref:`example_plot_rfe_digits.py`,
             :ref:`example_plot_feature_selection.py`

"""
print(__doc__)

import numpy as np
from sklearn.feature_selection import RFE, SelectKBest, chi2, \
    MinRedundancyMaxRelevance
from sklearn.linear_model import LogisticRegression


#  Number of samples in the dataset
N = 100

#  Associating a class to each sample in the dataset
y = np.array([0] * 50 + [1] * 50)

#  Creating a feature able to classify 60% of the samples
A = np.array([0] * 30 + [1] * 20 + [1] * 20 + [2] * 30)

#  Creating a feature able to classify 50% of the samples
B = np.array([2] * 25 + [1] * 25 + [1] * 25 + [0] * 25)

#  Creating a feature able to classify 40% of the samples
C = np.array([2] * 20 + [0] * 30 + [1] * 30 + [2] * 20)

X = np.array([A, B, C]).T
feature = ['A', 'B', 'C']

# We will be using the following three selectors
selectors = [('RFE', RFE(LogisticRegression(), 2)),
             ('Uni', SelectKBest(chi2, k=2)),
             ('mRMR', MinRedundancyMaxRelevance(k=2))]

for name, selector in selectors:
    k = selector.fit(X, y).get_support(True).tolist()
    print name, 'selected %s and %s' % (feature[k[0]], feature[k[1]])
