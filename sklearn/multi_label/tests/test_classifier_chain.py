from sklearn.multi_label import ClassifierChain
from sklearn.datasets import make_multilabel_classification
from sklearn.linear_model import LogisticRegression
from sklearn.utils.testing import assert_equal


def test_classifier_chain_fit_and_predict():
    """fit a model and assert that the label predictions have the same shape as the label targets"""

    X, Y = make_multilabel_classification(n_samples=10000, n_features=100, n_classes=10)
    classifier_chain = ClassifierChain(LogisticRegression())
    classifier_chain.fit(X, Y)
    Y_pred = classifier_chain.predict(X)
    assert_equal(Y_pred.shape, Y.shape)
