"""
The :mod:`sklearn.metrics.cluster` submodule contains evaluation metrics for
cluster analysis results. There are two forms of evaluation:

- supervised, which uses a ground truth class values for each sample.
- unsupervised, which does not and measures the 'quality' of the model itself.
"""
from .supervised import adjusted_mutual_info_score
from .supervised import adjusted_rand_score
from .supervised import completeness_score
from .supervised import contingency_matrix
from .supervised import expected_mutual_information
from .supervised import homogeneity_completeness_v_measure
from .supervised import homogeneity_score
from .supervised import mutual_info_score
from .supervised import v_measure_score
from .unsupervised import silhouette_samples
from .unsupervised import silhouette_score
