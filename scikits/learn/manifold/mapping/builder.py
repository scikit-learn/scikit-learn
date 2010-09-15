# -*- coding: utf-8 -*-

import inspect

from ...base import BaseEstimator

from barycenter import Barycenter

def builder(embedding, kind = "Barycenter", n_neighbors = None, neigh = None,
    neigh_alternate_arguments = None):
    """
    Function that will create a builder depending on the arguments it is passed
    
    Parameters
    ----------
    embedding :
        A usable embedding instance
    
    kind : object
        The type of mapper to use. Can be:
            * None : no mapping built
            * "Barycenter" (default) : Barycenter mapping
            * a class object : a class that will be instantiated with the
                arguments of this function
            * an instance : an instance that will be fit() and then
                transform()ed

    n_neighbors : int
      The number of K-neighboors to use (optional, default 9) if neigh is not
      given.

    neigh : Neighbors
      A neighboorer (optional). By default, a K-Neighbor research is done.
      If provided, neigh must be a functor class . `neigh_alternate_arguments` 
      will be passed to this class constructor.

    neigh_alternate_arguments : dictionary
      Dictionary of arguments that will be passed to the `neigh` constructor
    """
    if kind == None:
        return None
    elif kind == "Barycenter":
        mapping = Barycenter(n_neighbors = n_neighbors, neigh = neigh,
        neigh_alternate_arguments = neigh_alternate_arguments)
    elif isinstance(kind, BaseEstimator):
        mapping = kind
    elif inspect.isclass(kind) and issubclass(kind, BaseEstimator):
        mapping = kind(n_neighbors = n_neighbors, neigh = neigh,
            neigh_alternate_arguments = neigh_alternate_arguments)
    else:
      raise RuntimeError(
          "Argument 'kind' should be a BaseEstimator subclass or instance")
    mapping.fit(embedding)
    return mapping
