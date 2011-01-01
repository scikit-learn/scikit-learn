# -*- coding: utf-8 -*-

from ...neighbors import Neighbors

class NewNeighbors(Neighbors):
  predict = Neighbors.kneighbors

  def __init__(self, k, **kwargs):
    Neighbors.__init__(self, k)
