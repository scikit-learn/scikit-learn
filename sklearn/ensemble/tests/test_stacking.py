"""
Testing for the stacking ensemble module (sklearn.ensemble.stacking).
"""

# Author: Caio Oliveira
# License BSD 3 clause

from sklearn.utils.testing import (assert_equal, assert_array_equal,
                                   assert_almost_equal)
from sklearn.ensemble import (BlendedClassifierTransformer, make_stack_layer,
                              make_stacked_classifier)
from sklearn.linear_model import (LogisticRegression, RidgeClassifier)
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn import datasets

iris = datasets.load_iris()
X, y = iris.data[:, 1:3], iris.target


def test_transformer_init():
    # Tests if attributes are set correctly

    base_clf = RidgeClassifier(random_state=1)
    blended = BlendedClassifierTransformer(base_clf)

    assert_equal(blended.clf, base_clf)
    assert_equal(blended.cv, 3)


def test_stacking_api():
    # Just test if the API generates the expected pipelines
    clf1 = RidgeClassifier(random_state=1)
    clf2 = RidgeClassifier(random_state=1)
    clf3 = RidgeClassifier(random_state=1)
    clf4 = RidgeClassifier(random_state=1)
    clf5 = RidgeClassifier(random_state=1)

    layer = make_stack_layer(clf1, clf2, clf3)

    assert_array_equal([x[1].clf for x in layer.transformer_list],
                       [clf1, clf2, clf3])

    clf_matrix = [[clf1, clf2], [clf3, clf4]]
    final_clf = make_stacked_classifier(clf_matrix, clf5)

    matrix_from_pipeline = [[x[1].clf for x in y[1].transformer_list]
                            for y in final_clf.steps[:-1]]

    assert_array_equal(clf_matrix, matrix_from_pipeline)

    assert_equal(final_clf.steps[-1][1], clf5)


def test_classification():
    # Tests classification
    clf = make_stacked_classifier([[RandomForestClassifier(random_state=1),
                                    RandomForestClassifier(random_state=2)],
                                   [RandomForestClassifier(random_state=3),
                                    RandomForestClassifier(random_state=4)]],
                                  LogisticRegression(random_state=5))

    clf.fit(X, y)

    assert_almost_equal(.97333333334, accuracy_score(y, clf.predict(X)))
