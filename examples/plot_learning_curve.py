import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import PassiveAggressiveClassifier
from sklearn.datasets import load_digits
from sklearn.learning_curve import learning_curve

if __name__ == "__main__":
    estimator = PassiveAggressiveClassifier(n_iter=1)
    digits = load_digits()
    X, y = digits.data, digits.target

    plt.title("Learning Curves (Passive-Aggressive Classifier on Digits)")
    plt.xlabel("Training examples")
    plt.ylabel("Score")

    n_samples_range, train_scores, test_scores = learning_curve(
        estimator, X, y, cv=10, exploit_incremental_learning=False,
        n_jobs=1, verbose=False)
    plt.plot(n_samples_range, train_scores, label="Training score")
    plt.plot(n_samples_range, test_scores, label="Cross-validation score")

    n_samples_range, train_scores, test_scores = learning_curve(
        estimator, X, y, cv=10, exploit_incremental_learning=True,
        n_jobs=1, verbose=False)
    plt.plot(n_samples_range, train_scores, label="Training score")
    plt.plot(n_samples_range, test_scores, label="Cross-validation score")

    plt.legend(loc="best")
    plt.show()
