At each decision node in the tree, candidate split points are evaluated for each feature. For numerical features, the candidate thresholds are determined by sorting the unique values observed in the training data and using the midpoints between consecutive unique values as potential split thresholds. If a feature has `n` unique values, then `n - 1` possible splits are considered. This exhaustive approach affects the computational cost per node, which grows linearly with the number of unique values per feature.

The features are always randomly permuted at each split
