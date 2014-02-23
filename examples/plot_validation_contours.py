import matplotlib.pyplot as plt
import numpy as np
from sklearn.datasets import load_digits
from sklearn.svm import SVC
from sklearn.learning_curve import validation_curve
from sklearn.externals.joblib import Memory

memory = Memory(cachedir=".", verbose=0)

@memory.cache
def grid(X, y, Cs, gammas):
    param_grid = {"C": Cs, "gamma": gammas}

    tr, te, times = validation_curve(SVC(kernel="rbf"), X, y, param_grid, cv=3)

    shape = (len(Cs), len(gammas))
    tr = tr.mean(axis=1).reshape(shape)
    te = te.mean(axis=1).reshape(shape)
    times = times.mean(axis=1).reshape(shape)

    return tr, te, times

digits = load_digits()
X, y = digits.data, digits.target

gammas = np.logspace(-6, -1, 5)
Cs = np.logspace(-3, 3, 5)

tr, te, times = grid(X, y, Cs, gammas)


for title, values in (("Training accuracy", tr),
                      ("Test accuracy", te),
                      ("Training time", times)):

    plt.figure()

    plt.title(title)
    plt.xlabel("C")
    plt.xscale("log")

    plt.ylabel("gamma")
    plt.yscale("log")

    X1, X2 = np.meshgrid(Cs, gammas)
    cs = plt.contour(X1, X2, values)

    plt.colorbar(cs)

plt.show()
