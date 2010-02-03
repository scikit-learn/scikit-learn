
# Matthieu Brucher
# Last Change : 2008-04-07 15:47

"""
PCA module, by Zachary Pincus
"""

import numpy
import numpy.linalg

def PCA(samples, nb_coords, **kwargs):
  """
  Performs a PCA data reduction
  """
  centered = samples - numpy.mean(samples, axis=0)
  try:
    corr = numpy.dot(centered.T, centered)
    (w, v) = numpy.linalg.eigh(corr)
    index = numpy.argsort(w)

    unscaled = v[index[-1:-1-nb_coords:-1]]
    vectors = unscaled#numpy.sqrt(w[index[-1:-1-nb_coords:-1]])[:,numpy.newaxis] * unscaled
    inv = numpy.linalg.inv(numpy.dot(unscaled.T, unscaled.T))
    return numpy.dot(centered, numpy.dot(vectors.T, inv))

  except:
    corr = numpy.dot(centered, centered.T)
    (w, v) = numpy.linalg.eigh(corr)
    index = numpy.argsort(w)

    unscaled = v[:,index[-1:-1-nb_coords:-1]]
    vectors = numpy.dot(unscaled.T, centered)
    vectors = (1/numpy.sqrt(w[index[-1:-1-nb_coords:-1]])[:,numpy.newaxis]) * vectors
    return numpy.dot(centered, vectors.T)