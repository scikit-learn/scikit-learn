import numpy as np
import matplotlib.pyplot as plt
from sklearn.naive_bayes import GaussianNB
from sklearn.datasets import load_digits
from sklearn.learning_curve import learning_curve # TODO should be: from sklearn import learning_curve

if __name__ == "__main__":
    estimator = GaussianNB()
    digits = load_digits()
    X, y = digits.data, digits.target

    n_samples_range, train_scores, test_scores = learning_curve(
        estimator, X, y, cv=10, n_jobs=4, verbose=False)

    plt.title("Learning Curves (Naive Bayes on Digits Dataset)")
    plt.xlabel("Training examples")
    plt.ylabel("Score")
    plt.plot(n_samples_range, train_scores, label="Training score")
    plt.plot(n_samples_range, test_scores, label="Cross-validation score")
    plt.legend(loc="best")
    plt.show()
