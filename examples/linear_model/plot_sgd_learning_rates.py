"""
==================================
Comparing learning rates of SGD
==================================

An example showing how different learning_rate
schedules affect accuracy convergence

"""
# Author: Danny Sullivan <dsullivan7 at hotmail dot com>

import numpy as np
import matplotlib.pyplot as plt
from sklearn import datasets
from sklearn.linear_model import SGDClassifier
from sklearn.cross_validation import train_test_split

epochs = 20
train_proportion = .8


def train(X_train, y_train, X_test, y_test, classes, classifiers):
    plt.figure()
    for name, clf, scores in classifiers:
        print("training", name, "...")
        for r in range(epochs):
            clf.partial_fit(X_train, y_train, classes=classes)
            score = clf.score(X_test, y_test)
            scores.append(1.0 - score)

        plt.plot(1 + np.arange(epochs), scores, label=name)
        print(name, scores)

    plt.legend(loc="upper right")
    plt.xlabel("epoch")
    plt.ylabel("error")


if __name__ == "__main__":
    alpha = 1.e-5
    news_clfs = [
        ("invscaling", SGDClassifier(learning_rate="invscaling",
                                     alpha=alpha,
                                     n_iter=20,
                                     eta0=50.0), []),
        ("optimal", SGDClassifier(learning_rate="optimal",
                                  alpha=alpha,
                                  n_iter=20), []),
        # ("adadelta", SGDClassifier(learning_rate="adadelta",
        #                            alpha=alpha,
        #                            n_iter=20,
        #                            eps0=0.2,
        #                            rho0=.9), []),
        ("adadelta", SGDClassifier(learning_rate="adadelta",
                                   alpha=alpha,
                                   n_iter=20,
                                   eps0=100.0,
                                   rho0=.8), []),
    ]

    alpha = 1.e-8
    digits_clfs = [
        ("invscaling", SGDClassifier(learning_rate="invscaling",
                                     alpha=alpha,
                                     n_iter=20,
                                     eta0=.01), []),
        ("optimal", SGDClassifier(learning_rate="optimal",
                                  alpha=alpha,
                                  n_iter=20), []),
        ("adagrad", SGDClassifier(learning_rate="adagrad",
                                  alpha=alpha,
                                  n_iter=20,
                                  eta0=.01,
                                  eps0=.0), []),
    ]

    alpha = 1.e-5
    faces_clfs = [
        # ("invscaling", SGDClassifier(learning_rate="invscaling",
        #                              alpha=alpha,
        #                              n_iter=20,
        #                              eta0=40.), []),
        # ("optimal", SGDClassifier(learning_rate="optimal",
        #                           alpha=alpha,
        #                           n_iter=20), []),
        # ("adagrad", SGDClassifier(learning_rate="adagrad",
        #                           alpha=alpha,
        #                           n_iter=20,
        #                           eta0=50.0,
        #                           eps0=0.1), []),
        # ("adadelta", SGDClassifier(learning_rate="adadelta",
        #                            alpha=alpha,
        #                            n_iter=20,
        #                            eps0=1.0,
        #                            rho0=.95), []),
    ]

    tests = [
        (datasets.load_digits(), digits_clfs),
        (datasets.fetch_20newsgroups_vectorized(), news_clfs),
        # (datasets.fetch_lfw_pairs(), faces_clfs),
    ]

    for data_set, classifiers in tests:
        X, y = data_set.data, data_set.target
        classes = np.unique(y)

        train_num = int(train_proportion * X.shape[0])

        X_train, X_test, y_train, y_test = \
            train_test_split(X, y, train_size=train_proportion,
                             random_state=77)

        train(X_train, y_train, X_test, y_test, classes, classifiers)

    plt.show()








