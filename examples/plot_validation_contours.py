import matplotlib.pyplot as plt
import numpy as np
from sklearn.datasets import load_digits
from sklearn.svm import SVC
from sklearn.learning_curve import validation_curve
from sklearn.externals.joblib import Memory

memory = Memory(cachedir=".", verbose=0)

@memory.cache
def grid(X, y, Cs, gammas):
    scores = np.zeros((len(Cs), len(gammas)))

    for i, C in enumerate(Cs):
        tr, te = validation_curve(SVC(kernel="rbf", C=C), X, y,
                                  param_name="gamma", param_range=gammas, cv=3)
        scores[i] = te.mean(axis=1)

    return scores

digits = load_digits()
X, y = digits.data, digits.target

gammas = np.logspace(-6, -1, 5)
Cs = np.logspace(-3, 3, 5)

scores = grid(X, y, Cs, gammas)


plt.xlabel("C")
plt.xscale("log")

plt.ylabel("gamma")
plt.yscale("log")

X1, X2 = np.meshgrid(Cs, gammas)
cs = plt.contour(X1, X2, scores)

plt.colorbar(cs)

plt.show()
