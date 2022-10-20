"""
For many to many matching!
For example, if you have n products from one source and m products from another
you can vectorize them with the same number of features and use this function for matching
"""


import numpy as np

def feature_matching(nd1, nd2):
    """:param nd1 is an array of shape n*x
        :param nd2 is an array of shape m*x
        :return an array of shape n*m all zero except for the matching rows
        if row i of nd1 and row j of nd2 are equal the value of the returned matrix at position [i,j] is 1"""
    assert nd1.shape[1] == nd2.shape[1]
    inds = []
    v1 = np.repeat(nd1, repeats=nd2.shape[0], axis=0)
    v2 = np.array([nd2] * nd1.shape[0]).reshape(nd1.shape[0] * nd2.shape[0], nd1.shape[1])
    ms = np.where(np.all(v1 == v2, axis=1))[0]
    for id in ms:
        inds.append((int(id / nd2.shape[0]), id % nd2.shape[0]))
    vz = np.zeros((nd1.shape[0], nd2.shape[0]))
    rows, cols = zip(*inds)
    vz[rows, cols] = [1] * len(inds)
    return vz