import numpy as np
from sklearn.tree import _test


def main():
    crit = _test.Criterion()
    splitter = _test.Splitter(crit, 0)

    X = np.array([[0., 1], [0, 0]])
    splitter.test(X)


if __name__ == '__main__':
    main()