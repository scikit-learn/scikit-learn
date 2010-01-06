
"""
Piecewise Linear Mapping Regression module
"""

# Matthieu Brucher
# Last Change : 2008-11-06 10:50

import math
import numpy
import numpy.linalg as linalg
from numpy.random import shuffle
from scikits.learn.machine.manifold_learning import stats

class PLMR(object):
  """
  Regression with piecewise linear functions
  Uses ML or mean square error (same error for every piecewise function
  """
  def __init__(self, points, coords, neighbors, random_variable = stats.IsotropicGaussianVariable, RBF_field = stats.RBFField, correction_factor = 7.0):
    """
    Initializes the regression
    - points are the initial points
    - coords are the coordinates that will be used
    - neighbors is the number of neighbour used for determining a plan's equation
    - random_variable is the kid of random variable that will be used for estimation, it is supposed to be identical for every piecewise function
    - correction_factor is the factor for belonging setting
    """
    self.points = points
    self.coords = numpy.append(coords, numpy.ones((len(coords),1)), axis = 1).copy()
    self.correction_factor = correction_factor
    self.graph = self.create_graph(coords, neighbors)

    self.random_variable = random_variable()
    self.RBF_field = RBF_field
    self.coords_field = RBF_field(coords, weight = 1)

  def learn(self):
    """
    Tries to learn the model
    """
    self.belonging_vector = numpy.ones(len(self.coords), dtype = numpy.int) * -1
    self.equations = []

    self.findEquations()
    self.findEquations(False)
    self.assignOutliers()

    self.computeError()
    del self.points
    self.RBFF = [self.createRBF(numpy.where(self.belonging_vector == plan)[0]) for plan in range(0, len(self.equations))]

  def __getstate__(self):
    return (self.coords, self.random_variable, self.RBF_field, self.equations, self.belonging_vector)

  def __setstate__(self, state):
    self.coords = state[0]
    self.random_variable = state[1]
    self.RBF_field = state[2]
    self.coords_field = self.RBF_field(numpy.array(self.coords[:,:-1]), weight = 1)
    self.equations = state[3]
    self.belonging_vector = state[4]
    self.RBFF = [self.createRBF(numpy.where(self.belonging_vector == plan)[0]) for plan in range(0, len(self.equations))]

  def create_graph(self, coords, neighbors):
    """
    Creates a pseudo graph of the nearest neighbors
    """
    import neighbors as tool_neighbors

    graph = [set() for i in xrange(len(coords))]
    self.neighbourer = tool_neighbors.Kneighbors(coords, neighbors)

    for point in range(0, len(coords)):
      for neighbour in self.neighbourer(coords[point]):
        graph[point].add(neighbour[1])
        graph[neighbour[1]].add(point)

    return [list(neighbors) for neighbors in graph]

  def findEquations(self, random = True):
    """
    Tries to find randomly equations in the space
    """
    order = numpy.arange(0, len(self.coords) - 1)
    if random:
      shuffle(order)
    for value in order:
      if (self.belonging_vector[self.graph[value]] == -1).all():
        self.findEquationAround(value)
        self.updateBV()
        self.pruneEquations()
        self.updateEquations()
        self.updateBV()
        self.pruneEquations()

  def findEquationAround(self, value):
    """
    Tries to find, if possible, a plan around the point indicated by value
    """
    self.equations.append(self.computeEquation2(self.coords[self.graph[value]], self.points[self.graph[value]]))
    self.belonging_vector[self.graph[value]] = len(self.equations) - 1

  def computeEquation2(self, coords, points):
    """
    Computes an equation from its coordinates and the corresponding points
    """
    coords_bis = numpy.asmatrix(coords)
    points_bis = numpy.asmatrix(points)
    equation = linalg.inv(coords_bis.T*coords_bis).T * coords_bis.T * points_bis
    return numpy.asarray(equation)

  def updateBV(self):
    """
    Updates the belonging vector
    """
    if self.equations == []:
      self.belonging_vector[:] = -1
      return

    errors = [numpy.sum((self.points - numpy.dot(self.coords, equation))**2, axis=1) for equation in self.equations]
    residuals = self.computeResiduals()
    variance = numpy.mean(numpy.sum(residuals**2, axis=1))

    errors = numpy.array(errors)
    best = errors.argmin(axis = 0)
    corr_best = errors.min(axis = 0)
    validated = numpy.where(corr_best < self.correction_factor * variance)
    self.belonging_vector[validated] = best[validated]

  def pruneEquations(self):
    """
    If an equation does not have enough points, it is deleted
    """
    min_size = (self.coords.shape[1] - 1) * 2
    for plan in range(max(self.belonging_vector), -1, -1):
      if len(numpy.where(self.belonging_vector == plan)[0]) < min_size:
        if plan < len(self.equations):
          del self.equations[plan]
        self.belonging_vector[numpy.where(self.belonging_vector == plan)[0]] = -1
        self.belonging_vector[numpy.where(self.belonging_vector > plan)[0]] -= 1

  def updateEquations(self):
    """
    Updates plan equations
    """
    size = numpy.max(self.belonging_vector)
    self.equations = []
    for plan in range(0, size + 1):
      self.equations.append(self.computeEquation2(self.coords[numpy.where(self.belonging_vector == plan)[0]], self.points[numpy.where(self.belonging_vector == plan)[0]]))

  def assignOutliers(self):
    """
    Assign remaining outliers
    """
    for outlier in numpy.where(self.belonging_vector == -1)[0]:
      self.belonging_vector[outlier] = self.computeNearestPlan(self.coords[outlier])

  def computeError(self):
    """
    Computes the mean and variance error
    """
    residuals = self.computeResiduals()
    self.random_variable.addSample(residuals)
    self.random_variable.compute()
    self.random_variable.clean()

  def computeNearestPlan(self, coords):
    """
    Returns the index of the nearest plan
    """
    RBFFs = [self.createRBF(numpy.where(self.belonging_vector == plan)[0]) for plan in range(0, len(self.equations))]
    p = [RBFF(coords[:-1]) for RBFF in RBFFs]
    return p.index(max(p))

  def computeResiduals(self):
    """
    Computes all residuals
    """
    if len(self.equations) > 0:
      errors = [(self.points - numpy.dot(self.coords, equation)) for equation in self.equations]
      residuals = numpy.zeros(self.points.shape)
      for plan in range(0, len(self.equations)):
        residuals[numpy.where(self.belonging_vector == plan)[0]] = errors[plan][numpy.where(self.belonging_vector == plan)[0]]
      return residuals[numpy.where(self.belonging_vector >= 0)[0]]
    else:
      return self.points

  def createRBF(self, indexes):
    """
    Creates an RBF around some points diven by their indexes
    """
    weight = float(len(indexes))/len(self.belonging_vector)
    return self.RBF_field(self.coords[indexes,:-1], weight = weight)

  def ensure_connexity(self):
    """
    Ensure that every set of points is connected, each set being the set of points labeled to a plan
    Return True if the equations must be computed again
    """
    nbPlans = len(self.equations)

    for plan in range(0, nbPlans):
      coords = set(numpy.where(self.belonging_vector == plan)[0])
      while coords != set():
        el = coords.pop()
        coords.add(el)
        component = self.find_component(el, coords)
        coords.difference_update(component)
        if coords != set():
          print coords
          self.belonging_vector[list(component)] = nbPlans
          nbPlans += 1
    if nbPlans != len(self.equations):
      print "Connexity made plans be split"
      return True
    return False

  def find_component(self, el, coords):
    """
    Find the component in the subgraph coords where el is
    """
    component = set([el])
    component2 = set()
    while component != component2:
      component2 = component
      component = set()
      for point in component2:
        component.update(self.graph[point])
      component.intersection_update(coords)

    return component

  def get_log_likelihood(self, coords, point, mask=1., **kwargs):
    """
    Returns the negative log-likelihood for a given point and set of coordinates
    """
    if 'equation' in kwargs:
      equation = kwargs['equation']
    else:
      equation = self.computeNearestPlan(coords)

    reconstruct = numpy.dot(coords, self.equations[equation])
    epsilon = point - reconstruct

    return - self.random_variable.RBF(epsilon)

  def get_MAP(self, coords, points, mask=1., **kwargs):
    """
    Returns the MAP for a given point
    """
    if 'equation' in kwargs:
      somme = - self.coords_field(coords[:-1])
    else:
      somme = - self.RBFF[self.computeNearestPlan(coords)](coords[:-1])
    cost = self.get_log_likelihood(coords, points, mask, **kwargs)
    return cost + somme

  def get_point(self, coords):
    """
    Computes a point based on its coordinates
    """
    equation = self.computeNearestPlan(coords)
    return numpy.dot(coords, self.equations[equation])

  def get_random_args(self):
    return self.random_variable.get()

  def set_random_args(self, args):
    return self.random_variable.set(args)

  def setup_random(self):
    return self.random_variable.setup()
