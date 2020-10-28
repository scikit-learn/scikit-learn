from numpy import ndarray


class VoronoiINode:
	def __init__(self, data_size: int, split_points: ndarray, children: ndarray, depth: int, node_index: int):
		self.data_size: int = data_size
		self.split_points: ndarray = split_points
		self.children: ndarray = children
		self.depth: int = depth
		self.node_index: int = node_index
