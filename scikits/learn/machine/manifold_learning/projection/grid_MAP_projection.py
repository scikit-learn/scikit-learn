
"""
Projection with MAP on a piecewise linear function module with a grid
"""

# Matthieu Brucher
# Last Change : 2008-03-13 16:08

from grid_ML_projection import *

__all__ = ['GridMAPProjection']

class GridMAPProjection(GridMLProjection):
  """
  Class that will handle the projection
  - PLMR is an instance of PLMR or that satisfies its attribute interface
  """
  def __init__(self, PLMR):
    GridMLProjection.__init__(self, PLMR)
    self.PLMRcost = self.PLMR.get_MAP
