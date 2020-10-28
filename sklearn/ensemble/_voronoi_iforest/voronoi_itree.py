from .voronoi_inode import VoronoiINode
from numpy import arange, dot, empty, euler_gamma, eye, flatnonzero, full, log, ndarray, split, sum, where, zeros_like
from numpy.random import choice
from scipy.spatial.distance import cdist


class VoronoiITree:
	def __init__(self, branching_factor: int, metric: str = 'euclidean'):
		self.branching_factor: int = branching_factor
		self.metric: str = metric
		self.depth_limit: float = 0.0
		self.root: VoronoiINode = None
		self.nodes_count: int = 0

	def build(self, data: ndarray) -> 'VoronoiITree':
		# Adjust depth limit according to data cardinality and branching factor
		self.depth_limit: float = log(data.shape[0]) / log(self.branching_factor)
		# Recursively build the tree
		self.nodes_count, self.root = self.recursive_build(data)
		return self

	def recursive_build(self, data: ndarray, depth: int = 0, node_index: int = 0) -> (int, VoronoiINode):
		# If there aren't enough samples to be split according to the branching factor or the depth limit has been
		# reached, build a leaf node
		if data.shape[0] < self.branching_factor or depth >= self.depth_limit:
			return node_index + 1, VoronoiINode(data.shape[0], None, None, depth, node_index)
		else:
			# Generate split points
			split_points: ndarray = data[choice(data.shape[0], size=self.branching_factor, replace=False)]
			# Partition data
			partition_indices: list = self.split_data(data, split_points)
			# Generate recursively children nodes
			children: ndarray = empty(shape=(self.branching_factor,), dtype=VoronoiINode)
			for i, indices in enumerate(partition_indices):
				node_index, children[i] = self.recursive_build(data[indices], depth + 1, node_index)
			return node_index + 1, VoronoiINode(data.shape[0], split_points, children, depth, node_index)

	def predict(self, data: ndarray) -> ndarray:
		# Compute depth of each sample
		return self.recursive_depth_search(self.root, data, empty(shape=(data.shape[0],), dtype=float))

	def recursive_depth_search(self, node: VoronoiINode, data: ndarray, depths: ndarray) -> ndarray:
		# If the current node is a leaf, fill the depths vector with the current depth plus a normalization factor
		if node.split_points is None:
			depths[:] = node.depth + self.get_random_path_length(node.data_size)
		else:
			# Partition data
			partition_indices: list = self.split_data(data, node.split_points)
			# Fill the vector of depths
			for i, indices in enumerate(partition_indices):
				depths[indices]: ndarray = self.recursive_depth_search(node.children[i], data[indices], depths[indices])
		return depths

	def apply(self, data: ndarray) -> ndarray:
		# Compute leaf index of each sample
		return self.recursive_index_search(self.root, data, empty(shape=(data.shape[0],), dtype=int))

	def recursive_index_search(self, node: VoronoiINode, data: ndarray, leaves_index: ndarray) -> ndarray:
		# If the current node is a leaf, fill the leaves index vector with the current node index
		if node.split_points is None:
			leaves_index[:] = node.node_index
		else:
			# Partition data
			partition_indices: list = self.split_data(data, node.split_points)
			# Fill the vector of leaves index
			for i, indices in enumerate(partition_indices):
				leaves_index[indices]: ndarray = self.recursive_index_search(node.children[i], data[indices], leaves_index[indices])
		return leaves_index

	def weight(self, data: ndarray) -> ndarray:
		# Compute leaf mass of each sample
		return self.recursive_mass_search(self.root, data, empty(shape=(data.shape[0],), dtype=int))

	def recursive_mass_search(self, node: VoronoiINode, data: ndarray, masses: ndarray) -> ndarray:
		# If the current node is a leaf, fill the masses vector with the current node mass
		if node.split_points is None:
			masses[:] = node.data_size
		else:
			# Partition data
			partition_indices: list = self.split_data(data, node.split_points)
			# Fill the vector of leaves index
			for i, indices in enumerate(partition_indices):
				masses[indices]: ndarray = self.recursive_mass_search(node.children[i], data[indices], masses[indices])
		return masses

	def decision_path(self, data: ndarray) -> ndarray:
		# Compute path of each sample
		return self.recursive_path_search(self.root, data, full((data.shape[0], self.nodes_count), False))

	def recursive_path_search(self, node: VoronoiINode, data: ndarray, paths: ndarray) -> ndarray:
		# Fill the position of the actual node
		paths[:, node.node_index] = True
		# If the current node is not a leaf,
		if node.split_points is not None:
			# Partition data
			partition_indices: list = self.split_data(data, node.split_points)
			# Fill the vector of leaves index
			for i, indices in enumerate(partition_indices):
				paths[indices]: ndarray = self.recursive_path_search(node.children[i], data[indices], paths[indices])
		return paths

	def split_data(self, data: ndarray, split_points: ndarray) -> list:
		# Compute distances of data from split points
		if self.metric == 'tanimoto':
			distances: ndarray = VoronoiITree.tanimoto_distance(data, split_points)
		else:
			distances: ndarray = cdist(data, split_points, metric=self.metric)
		# Build full membership mask
		full_membership: ndarray = distances == distances.min(axis=1, keepdims=True)
		# Keep only the first membership for each sample
		membership: ndarray = zeros_like(full_membership)
		membership[arange(full_membership.shape[0]), full_membership.argmax(axis=1)]: ndarray = True
		# Split data according to their membership
		row, col = where(membership.T)
		partition: list = split(col, flatnonzero(row[1:] != row[:-1]) + 1)
		return partition

	def get_random_path_length(self, num_samples: int) -> float:
		if self.branching_factor == 2:
			if num_samples <= 1:
				return 0
			elif num_samples == 2:
				return 1
			else:
				return 2.0 * (log(num_samples - 1.0) + euler_gamma) - 2.0 * (num_samples - 1.0) / num_samples
		else:
			if num_samples == 1:
				return 0
			elif 1 < num_samples <= self.branching_factor:
				return 1
			else:
				return log(num_samples) / log(self.branching_factor)

		# if num_samples == 1:
		# 	return 0
		# elif 1 < num_samples <= self.branching_factor:
		# 	return 1
		# else:
		# 	return (log(num_samples) + log(self.branching_factor - 1) + euler_gamma) / log(self.branching_factor) - 0.5

		# if num_samples < self.branching_factor:
		# 	return 0
		# elif num_samples == self.branching_factor:
		# 	return 1
		# else:
		# 	return log(num_samples) / log(self.branching_factor)

	@staticmethod
	def tanimoto_distance(XA: ndarray, XB: ndarray) -> ndarray:
		# Compute Tanimoto similarity for each couple of samples
		AB: ndarray = dot(XA, XB.T)
		A_squared: ndarray = sum(XA * XA, axis=1)
		B_squared: ndarray = sum(XB * XB, axis=1)
		Y: ndarray = AB / (A_squared[:, None] + B_squared[None, :] - AB)
		# Transform similarities into distances
		Y: ndarray = 1 - Y
		return Y
