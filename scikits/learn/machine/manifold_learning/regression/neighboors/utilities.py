"""
Neighboors utility functions
"""

# Matthieu Brucher
# Last Change : 2007-10-25 14:05

def create_graph(coords, neighboors):
  """
  Creates a pseudo graph of the nearest neighboors
  """
  import neighboors as neigh

  graph = [set() for i in xrange(len(coords))]
  neighboorer = neigh.KNeighboors(coords, neighboors)

  for point in range(0, len(coords)):
    for neighboor in neighboorer(coords[point]):
      graph[point].add(neighboor[1])
      graph[neighboor[1]].add(point)

  return [list(neighboors) for neighboors in graph]
