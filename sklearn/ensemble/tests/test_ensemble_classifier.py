"""Testing for the boost module (sklearn.ensemble.boost)."""

import numpy as np
from sklearn.utils.testing import assert_array_equal, assert_array_less
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises, assert_raises_regexp

from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB 
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import EnsembleClassifier
from sklearn import datasets

# Load the iris dataset and randomly permute it
iris = datasets.load_iris()
X, y = iris.data[:, 1:3], iris.target


def test_majority_label_iris():
    """Check classification by majority label on dataset iris."""
    
    np.random.seed(123)    
    clf1 = LogisticRegression()
    clf2 = RandomForestClassifier()
    clf3 = GaussianNB()
    eclf = EnsembleClassifier(clfs=[clf1, clf2, clf3])
    scores = cross_validation.cross_val_score(eclf, X, y, cv=5, scoring='accuracy')
    assert_almost_equal(scores.mean() == 0.94)
    

def test_weights_iris():
    """Check classification by average probabilities on dataset iris."""
    
    np.random.seed(123)    
    clf1 = LogisticRegression()
    clf2 = RandomForestClassifier()
    clf3 = GaussianNB()
    eclf = EnsembleClassifier(clfs=[clf1, clf2, clf3], weights=[1,2,10])
    scores = cross_validation.cross_val_score(eclf, X, y, cv=5, scoring='accuracy')
    assert_almost_equal(scores.mean() == 0.93)   




if __name__ == "__main__":
    import nose
    nose.runmodule()
