"""
:mod:`sklearn.metrics.cluster` is a module containing evaluation metrics for
cluster analysis results. There are two forms of evaluation:

- supervised, which uses a ground truth class values for each sample.
- unsupervised, which does not and measures the 'quality' of the model itself.
"""
from supervised import (homogeneity_completeness_v_measure,
                        homogeneity_score, completeness_score,
                        v_measure_score, adjusted_rand_score)
from unsupervised import silhouette_score, silhouette_samples
