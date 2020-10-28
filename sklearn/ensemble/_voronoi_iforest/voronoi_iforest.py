from .voronoi_itree import VoronoiITree
from concurrent.futures import (ThreadPoolExecutor,
                                wait)
from multiprocessing import cpu_count
from numpy import empty, euler_gamma, log, ndarray
from numpy.random import choice


class VoronoiIForest:
	def __init__(self, num_trees: int = 100, max_samples: int = 256, branching_factor: int = 2,
				 metric: str = 'euclidean', n_jobs: int = 1):
		self.num_trees: int = num_trees
		self.max_samples: int = max_samples
		self.branching_factor: int = branching_factor
		self.normalization_factor: float = None
		self.metric: str = metric
		self.trees: list = []
		self.n_jobs: int = cpu_count() if n_jobs == -1 else min(n_jobs, cpu_count())

	def fit(self, data: ndarray) -> None:
		# Clean the tree list
		self.trees: list = []
		# Adjust the number of samples to be picked according to data cardinality
		self.max_samples: int = min(self.max_samples, data.shape[0])
		# Adjust the branching factor according to max samples
		self.branching_factor: int = min(self.branching_factor, self.max_samples)
		# Compute the normalization factor
		self.normalization_factor: float = self.get_random_path_length(self.max_samples)
		# Build LSH trees
		with ThreadPoolExecutor(max_workers=self.n_jobs) as executor:
			futures: list = []
			for i in range(self.num_trees):
				sampled_data: ndarray = data[choice(data.shape[0], size=self.max_samples, replace=False)]
				futures.append(executor.submit(VoronoiITree(self.branching_factor, self.metric).build, sampled_data))
			wait(futures)
			self.trees: list = [future.result() for future in futures]

	def score_samples(self, data: ndarray) -> ndarray:
		# Compute the depths of all samples in each tree
		depths: ndarray = empty(shape=(data.shape[0], self.num_trees), dtype=float)
		with ThreadPoolExecutor(max_workers=self.n_jobs) as executor:
			futures: list = []
			for i, tree in enumerate(self.trees):
				futures.append(executor.submit(tree.predict, data))
			wait(futures)
			for i, future in enumerate(futures):
				depths[:, i] = future.result()
		# Compute the mean depth of each sample along all trees
		mean_depths: ndarray = depths.mean(axis=1)
		# Compute normalized mean depths
		normalized_mean_depths: ndarray = 2**(-mean_depths/self.normalization_factor)
		return normalized_mean_depths

	def apply(self, data: ndarray) -> ndarray:
		# Compute the leaves index of all samples in each tree
		indices: ndarray = empty(shape=(data.shape[0], self.num_trees), dtype=int)
		with ThreadPoolExecutor(max_workers=self.n_jobs) as executor:
			futures: list = []
			for i, tree in enumerate(self.trees):
				futures.append(executor.submit(tree.apply, data))
			wait(futures)
			for i, future in enumerate(futures):
				indices[:, i] = future.result()
		return indices

	def weight_samples(self, data: ndarray) -> ndarray:
		# Compute the weight of all samples in each tree
		weights: ndarray = empty(shape=(data.shape[0], self.num_trees), dtype=int)
		with ThreadPoolExecutor(max_workers=self.n_jobs) as executor:
			futures: list = []
			for i, tree in enumerate(self.trees):
				futures.append(executor.submit(tree.weight, data))
			wait(futures)
			for i, future in enumerate(futures):
				weights[:, i] = future.result()
		# Compute the mean weight of each sample along all trees
		mean_weights: ndarray = weights.mean(axis=1)
		return mean_weights

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
