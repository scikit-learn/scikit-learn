# -*- coding: utf-8 -*-

# Author: Christoph Hermes <hermes@hausmilbe.net>
#
# License: BSD Style.

import numpy as np
import scipy as sc

class PC:
    """
    Polynomial Classifier.
    """

    def __init__(self, degree=2):
        self.degree = degree
        pass

    def _build_pc_features(self, X):
        N,D = X.shape

        ## create polynomials up to desired degree
        numP = lambda n,p: sc.factorial(n+p)/(sc.factorial(n)*sc.factorial(p)) # number of polynomials
        # build up index list to combine features, -1 indicates unused feature 
        I = np.zeros((numP(D, self.degree), self.degree), dtype=int)-1
        for i in range(1, I.shape[0]):
            I[i,:] = I[i-1,:]
            for j in range(self.degree):
                if I[i-1,j]+1 < D:
                    I[i,j] = I[i-1,j]+1
                    break
            j -= 1
            while j>=0:
                I[i,j] = I[i,j+1]
                j -= 1

        # use index list to build combined polynomial features P
        P = np.ones((N, numP(D, self.degree)))
        for i in range(I.shape[0]):
            for d in range(self.degree):
                if I[i,d] > -1:
                    P[:,i] = P[:,i] * X[:,I[i,d]]

        return P

    def fit(self, X, y):
        N,D = X.shape

        self._classes = np.unique(y)

        ## create discriminant vector Y from y
        Y = np.zeros((N, len(self._classes)))
        for i in range(N):
            Y[i, self._classes==y[i]] = 1
        Y = np.matrix(Y)

        PX = np.matrix(self._build_pc_features(X))

        ## TODO: replace by Gauss-Jordan algorithm
        # direct computation of A
        self.A = np.linalg.inv(PX.T * PX / N) * (PX.T * Y / N)

        return self

    def predict(self, X):
        # build up polynomial features
        PX = np.matrix(self._build_pc_features(X))

        # apply classifier
        D = self.A.T * PX.T

        return np.array(np.argmax(D, axis=0))[0]

