"""
A comparison of different optimizers for NCA

Data: UC Irvine's Wine dataset.

"""
import numpy as np

from sklearn.base import clone
from sklearn import metrics
from sklearn.metric_learning import NCATransformer
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import Pipeline
from sklearn.cross_validation import KFold

from sklearn.utils.bench import total_seconds
from sklearn.datasets import fetch_mldata


def evaluate(est, X, y, scorer, k=10, train_size=0.7):
    kf = KFold(len(X), n_folds=k)
    scores = []
    
    for train_index, test_index in kf:
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        est = clone(est).fit(X_train, y_train)
        scores.append(scorer(est, X_test, y_test))
        
    return np.mean(scores)

if __name__ == '__main__':

    separator = "=" * 30
    wine_ds = fetch_mldata('wine')
    print("Loaded UCI Wine dataset")
    print(separator)

    knn1_nca = Pipeline([
                ('nca', NCATransformer(learning_rate=1, loss='kl', max_iter=100, solver='adagrad',
                                       random_state=0, verbose=2, n_init=3, method='semivectorized')),
                ('knn', KNeighborsClassifier(n_neighbors=1))
        ])

    print("Running 10-fold cross validation on Pipeline(NCA [semivectorized], 1NN)")
    nca_score = evaluate(knn1_nca, wine_ds.data, wine_ds.target, metrics.make_scorer(metrics.accuracy_score))
    print("Semivectorized NCA-assisted 1NN: accuracy_score = {:.5f}".format(nca_score))

    print(separator)

    knn1_nca = Pipeline([
                ('nca', NCATransformer(learning_rate=1, loss='kl', max_iter=100, solver='adagrad',
                                       random_state=0, verbose=2, n_init=3, method='vectorized')),
                ('knn', KNeighborsClassifier(n_neighbors=1))
        ])

    print("Running 10-fold cross validation on Pipeline(NCA [vectorized], 1NN)")
    nca_score = evaluate(knn1_nca, wine_ds.data, wine_ds.target, metrics.make_scorer(metrics.accuracy_score))
    print("Vectorized NCA-assisted 1NN: accuracy_score = {:.5f}".format(nca_score))

    print(separator)

    knn1 = KNeighborsClassifier(n_neighbors=1)
    print("Running 10-fold cross validation on 1NN")
    nca_score = evaluate(knn1, wine_ds.data, wine_ds.target, metrics.make_scorer(metrics.accuracy_score))
    print("1NN: accuracy_score = {:.5f}".format(nca_score))

