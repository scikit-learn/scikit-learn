
"""
Clustering regression
"""

# Matthieu Brucher
# Last Change : 2008-06-11 09:24

import numpy
import cluster

from MLPLMR import MLPLMR

class ClusteredRandomVariable(object):
  """
  A clustered random variable
  """
  def __init__(self, kind):
    """
    Initializes the clustered random variable
    """
    self.kind = kind
    self.PLMR = []

  def get(self):
    """
    Gets the parameters for each PLMR
    """
    l = []

    for PLMR in self.PLMR:
      l.append(PLMR.random_variable.get())
    return l

  def modify_kind(self, kind):
    """
    Modifies the RV used in the process
    """
    self.kind = kind
    for PLMR in self.PLMR:
      PLMR.random_variable = kind

  def set(self, args):
    """
    Sets the different random variables
    """
    for PLMR, arg in zip(self.PLMR, args):
      PLMR.random_variable.set(arg)

  def setup(self):
    """
    Sets up the different random variables
    """
    for PLMR in self.PLMR:
      PLMR.random_variable.setup()

class CPLMR(object):
  """
  Allows to compute a regression block-wise
  """
  def __init__(self, points, coords, neighbors, random_variable, RBF_field):
    """
    Initializes the regression
    - points are the initial points
    - coords are the coordinates that will be used
    - neighbors is the number of neighboor used for determining a plan's equation
    - random_variable is the kid of random variable that will be used for estimation, it is supposed to be identical for every piecewise function
    """
    self.points = points
    self.coords_orig = coords
    self.coords = numpy.append(coords, numpy.ones((len(coords),1)), axis = 1).copy()
    self.neighbors = neighbors

    self.RV = ClusteredRandomVariable(random_variable)

    self.random_variable = property(self.RV.modify_kind)

    self.RBF_field = RBF_field
    self.coords_field = RBF_field(coords, weight = 1)

  def learn(self):
    """
    Tries to learn the model
    """
    if not hasattr(self, 'clusters'):
      clusterer = cluster.GeneralCluster(0.33, 0.66)
      self.clusters = numpy.asarray(clusterer.process(numpy.corrcoef(self.points.T)))
      self.nbClusters = numpy.max(self.clusters)+1

    self.PLMR = []

    for cluster in range(0, self.nbClusters):
      points = self.points[:,numpy.where(self.clusters==cluster)[0]].copy()
      plmr = MLPLMR(points, self.coords_orig, self.neighbors, random_variable = self.RV.kind, RBF_field = self.RBF_field)
      plmr.learn()
      self.PLMR.append(plmr)

    self.RV.PLMR = self.PLMR

  def __getstate__(self):
    """
    Returns the state of the regression
    """
    return (self.RV, self.PLMR, self.clusters, self.coords_orig, self.coords_field)

  def __setstate__(self, state):
    """
    Sets the state of the regression
    """
    self.RV = state[0]
    self.PLMR = state[1]
    self.clusters = state[2]
    self.nbClusters = numpy.max(self.clusters)+1
    self.coords_orig = state[3]
    self.coords = numpy.append(self.coords_orig, numpy.ones((len(self.coords_orig),1)), axis = 1).copy()
    self.coords_field = state[4]
    self.random_variable = property(self.RV.modify_kind)

  def get_log_likelihood(self, coords, points, mask=1., **kwargs):
    """
    Returns the negative log-likelihood for a given point and set of coordinates
    """
    cost = 0

    for cluster in range(0, self.nbClusters):
      equation = self.PLMR[cluster].computeNearestPlan(coords)
      reconstruct = numpy.dot(coords, self.PLMR[cluster].equations[equation])
      epsilon = points[:,numpy.where(self.clusters==cluster)[0]] - reconstruct
      cost += self.PLMR[cluster].random_variable.RBF(epsilon)

    return - cost

  def get_MAP(self, coords, points, mask=1., **kwargs):
    """
    Returns the MAP for a given point
    """
    cost = self.get_log_likelihood(coords, points, mask, **kwargs)
    somme = - self.coords_field(coords[:-1])
    #somme = - sum([RBF(coords) for RBF in self.RBFF])
    return cost + somme

  def get_Point(self, coords):
    """
    Computes a point based on its coordinates
    """
    reconstruct = numpy.zeros(len(self.clusters))
    for cluster in range(0, self.nbClusters):
      equation = self.PLMR[cluster].computeNearestPlan(coords)
      reconstruct[numpy.where(self.clusters==cluster)[0]] = numpy.dot(coords, self.PLMR[cluster].equations[equation])

    return reconstruct
