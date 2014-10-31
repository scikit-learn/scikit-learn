"""
==================================
Comparing learning rates of SGD
==================================

An example showing how different learning_rate
schedules affect accuracy convergence

"""
# Author: Danny Sullivan <dsullivan7 at hotmail dot com>
#         Alex Gramfort

import numpy as np
import time
import matplotlib.pyplot as plt
from sklearn import datasets
from sklearn.linear_model import SGDClassifier
from sklearn.cross_validation import train_test_split

epochs = 20
train_proportion = .8


def train(X_train, y_train, X_test, y_test, classes, classifiers):
    plt.figure()
    colors = ['b', 'r', 'g', 'k', 'm']
    xx = 1 + np.arange(epochs)
    for col, (name, clf) in zip(colors, classifiers):
        scores_train = []
        scores_test = []
        print("Training %s ..." % name)
        t0 = time.time()
        for r in range(epochs):
            clf.partial_fit(X_train, y_train, classes=classes)
            scores_train.append(1.0 - clf.score(X_train, y_train))
            scores_test.append(1.0 - clf.score(X_test, y_test))

        print('Time elapsed: %s' % (time.time() - t0))
        plt.plot(xx, scores_train, col, label=name + ' (train)')
        plt.plot(xx, scores_test, col + '--', label=name + ' (test)')

    plt.legend(loc="upper right")
    plt.xlabel("epoch")
    plt.ylabel("error")
    plt.show()


if __name__ == "__main__":
    train_proportion = .8
    # shuffle = True
    shuffle = False
    loss = 'log'
    # loss = 'hinge'

    alpha = 1e-8
    news_clfs = [
        ("invscaling", SGDClassifier(learning_rate="invscaling",
                                     alpha=alpha, loss=loss,
                                     n_iter=20, shuffle=shuffle,
                                     eta0=60.0)),
        ("optimal", SGDClassifier(learning_rate="optimal",
                                  alpha=alpha, loss=loss,
                                  n_iter=20)),
        ("adadelta", SGDClassifier(learning_rate="adadelta",
                                   alpha=alpha, loss=loss,
                                   n_iter=20, shuffle=shuffle,
                                   eps0=100.,
                                   rho0=.8)),
        ("adagrad", SGDClassifier(learning_rate="adagrad",
                                  alpha=alpha, loss=loss,
                                  n_iter=20, shuffle=shuffle,
                                  eta0=60.,
                                  eps0=0.0)),
    ]

    alpha = 1e-1
    digits_clfs = [
        ("invscaling", SGDClassifier(learning_rate="invscaling",
                                     alpha=alpha, loss=loss,
                                     n_iter=20, shuffle=shuffle,
                                     eta0=.01)),
        ("optimal", SGDClassifier(learning_rate="optimal",
                                  alpha=alpha, loss=loss,
                                  n_iter=20)),
        ("adadelta", SGDClassifier(learning_rate="adadelta",
                                   alpha=alpha, loss=loss,
                                   n_iter=20, shuffle=shuffle,
                                   eps0=10.,
                                   rho0=.9)),
        ("adagrad", SGDClassifier(learning_rate="adagrad",
                                  alpha=alpha, loss=loss,
                                  n_iter=20, shuffle=shuffle,
                                  eta0=1.0,
                                  eps0=1e-2,
                                  fit_intercept=False)),
    ]

    alpha = 1e-1
    faces_clfs = [
        ("invscaling", SGDClassifier(learning_rate="invscaling",
                                     alpha=alpha, loss=loss,
                                     n_iter=20, shuffle=shuffle,
                                     eta0=1.)),
        ("optimal", SGDClassifier(learning_rate="optimal",
                                  alpha=alpha, loss=loss,
                                  n_iter=20)),
        ("adagrad", SGDClassifier(learning_rate="adagrad",
                                  alpha=alpha, loss=loss,
                                  n_iter=20, shuffle=shuffle,
                                  eta0=1.0,
                                  eps0=0.1)),
        ("adadelta", SGDClassifier(learning_rate="adadelta",
                                   alpha=alpha, loss=loss,
                                   n_iter=20, shuffle=shuffle,
                                   eps0=1.0,
                                   rho0=.95)),
    ]

    all_datasets = [
        (datasets.load_digits(), digits_clfs),
        # (datasets.fetch_20newsgroups_vectorized(), news_clfs),
        # (datasets.fetch_lfw_pairs(), faces_clfs),
    ]

    for data_set, classifiers in all_datasets:
        X, y = data_set.data, data_set.target
        classes = np.unique(y)

        train_num = int(train_proportion * X.shape[0])

        X_train, X_test, y_train, y_test = \
            train_test_split(X, y, train_size=train_proportion,
                             random_state=77)

        train(X_train, y_train, X_test, y_test, classes, classifiers)

    plt.show()








