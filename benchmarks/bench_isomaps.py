"""
=========================================
 Benchmark of Isomap and L-ISOMAP methods
=========================================

We will compare performance and evaluate the reductions of
example data sets using 1-NN classifier.
"""

from time import time

import numpy as np

from sklearn import datasets
from sklearn import manifold, neighbors
from sklearn.metrics import accuracy_score, r2_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline

n_samples = 1000
n_components = 2
landmarks = 'auto'
manifold_noise = .8
test_size = .2
n_jobs = 1
random_state = 0

classification_data_sets = (
    ('iris', datasets.load_iris()),
    ('digits', datasets.load_digits()),
    ('cancer', datasets.load_breast_cancer()),
)

regression_data_sets = (
    ('boston', datasets.load_boston()),
    ('diabetes', datasets.load_diabetes()),
)

manifold_data_sets = (
    ('s-curve', datasets.make_s_curve(n_samples)),
    ('swiss-roll', datasets.make_swiss_roll(n_samples)),
    ('noisy-s-curve', datasets.make_s_curve(n_samples, manifold_noise,
                                            random_state=random_state)),
    ('noisy-swiss-roll', datasets.make_swiss_roll(n_samples, manifold_noise,
                                                  random_state=random_state)),
)


def train_and_test_many(learner, data_sets, score_function):
    for ds_name, data_set in data_sets:
        print('==================')
        print('%s data set\n' % ds_name)

        isomap = manifold.Isomap(n_components=n_components)
        l_isomap = manifold.Isomap(n_components=n_components,
                                   landmarks=landmarks)

        i_pipe = Pipeline([('isomap', isomap), ('knn', learner)])
        l_pipe = Pipeline([('l_isomap', l_isomap), ('knn', learner)])

        for p_name, pipe in (
                ('Raw', learner),
                ('Isomap', i_pipe),
                ('L-Isomap', l_pipe)):
            if isinstance(data_set, tuple):
                data, target = data_set
            else:
                data, target = data_set.data, data_set.target

            X_train, X_test, y_train, y_test = train_test_split(
                data, target, test_size=test_size, random_state=0)

            train_dt = time()
            pipe.fit(X_train, y_train)
            train_dt = time() - train_dt

            y_pred = pipe.predict(X_test)

            print('%s (%.4f sec): %.4f'
                  % (p_name, train_dt, score_function(y_test, y_pred)))

        print()


if __name__ == '__main__':
    print(__doc__)

    np.random.seed(7434)

    print('Classification:')
    knn = neighbors.KNeighborsClassifier(n_jobs=n_jobs, n_neighbors=1)
    train_and_test_many(knn, classification_data_sets, accuracy_score)

    print('Regression:')
    knn = neighbors.KNeighborsRegressor(n_jobs=n_jobs, n_neighbors=1)
    train_and_test_many(knn, regression_data_sets, r2_score)

    print('Manifolds:')
    knn = neighbors.KNeighborsRegressor(n_jobs=n_jobs, n_neighbors=1)
    train_and_test_many(knn, manifold_data_sets, r2_score)
