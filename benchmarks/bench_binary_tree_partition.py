from sklearn.neighbors import KDTree
import numpy as np
from time import time


def main():
    test_cases = [
        ("ordered", np.arange(50000, dtype=float)),
        ("reverse ordered", np.arange(50000, 0, -1, dtype=float)),
        ("duplicated", np.zeros([50000], dtype=float))
    ]
    for name, case in test_cases:
        expanded_case = np.expand_dims(case, -1)
        begin = time()
        tree = KDTree(expanded_case, leaf_size=1)
        end = time()
        del tree
        print("{name}: {time}s".format(name=name, time=end - begin))


if __name__ == "__main__":
    main()
