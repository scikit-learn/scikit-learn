
"""
PCA Regression module
"""

# Matthieu Brucher
# Last Change : 2007-07-10 15:13

import math
import numpy
import numpy.linalg as linalg
from numpy.random import shuffle
import PLMR

class PCA(PLMR.PLMR):
  """
  Regression with one linear funtion
  """
  def learn(self):
    """
    Tries to learn the model
    """
    equation = numpy.zeros((self.coords.shape[1], self.points.shape[1]))
    equation[-1] = numpy.mean(self.points, axis=0)
    centered = self.points - equation[-1]

    coords_bis = numpy.asmatrix(self.coords[:,:-1])
    points_bis = numpy.asmatrix(centered)
    equation[0:-1] = linalg.inv(coords_bis.T*coords_bis).T * coords_bis.T * points_bis


    self.belonging_vector = numpy.zeros(len(self.coords), dtype = numpy.int)
    self.equations=[equation]

    self.computeError()
    del self.points
    self.RBFF = [self.createRBF(numpy.where(self.belonging_vector == plan)[0]) for plan in range(0, len(self.equations))]
