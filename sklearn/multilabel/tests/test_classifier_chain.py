from sklearn.multilabel import ClassifierChain
from sklearn.datasets import make_multilabel_classification
from sklearn.linear_model import LogisticRegression
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_not_equal
from sklearn.metrics import jaccard_similarity_score
from scipy import sparse


def test_fit_and_predict():
    """Fit classifier chain and verify performance"""
    X, Y = make_multilabel_classification(n_samples=10000,
                                          n_features=100,
                                          n_classes=10)
    classifier_chain = ClassifierChain(LogisticRegression())
    classifier_chain.fit(X, Y)
    Y_pred = classifier_chain.predict(X)
    Y_pred_binary = (Y_pred >= .5)
    assert_equal(Y_pred.shape, Y.shape)
    assert_greater(jaccard_similarity_score(Y, Y_pred_binary), 0.5)

def test_fit_and_predict_with_data_and_labels():
    """Fit classifier chain with sparse data and labels"""
    X, Y = make_multilabel_classification(n_samples=10000,
                                          n_features=100,
                                          n_classes=10)

    Y_sparse = sparse.csr_matrix(Y)
    X_sparse = sparse.csr_matrix(X)
    classifier_chain = ClassifierChain(LogisticRegression())
    classifier_chain.fit(X_sparse, Y_sparse)
    Y_pred = classifier_chain.predict(X_sparse)
    assert_equal(Y_pred.shape, Y.shape)


def test_order_shuffle():
    """Fit classifier chain with shuffled order"""
    X, Y = make_multilabel_classification(n_samples=10000,
                                          n_features=100,
                                          n_classes=10)
    classifier_chain = ClassifierChain(LogisticRegression(),
                                       order='random')
    classifier_chain.fit(X, Y)
    assert_not_equal(classifier_chain.order, list(range(10)))
    assert (len(classifier_chain.order) == 10)
    assert (len(set(classifier_chain.order)) == 10)


def test_classifiers_coef_size():
    """Fit classifier chain and verify the shape of coefficients"""
    X, Y = make_multilabel_classification(n_samples=10000,
                                          n_features=100,
                                          n_classes=10)
    classifier_chain = ClassifierChain(LogisticRegression(), shuffle=True)
    classifier_chain.fit(X, Y)

    assert_equal([c.coef_.size for c in classifier_chain.estimators_],
                 list(range(X.shape[1], X.shape[1] + Y.shape[1])))
