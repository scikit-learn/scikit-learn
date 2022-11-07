from collections import defaultdict
from itertools import product
from time import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.metrics import adjusted_mutual_info_score

repeat = 10
n_samples = [1_000, 10_000, 100_000, 1_000_000]
n_labels = [10, 100, 1_000]

rng = np.random.default_rng(0)

result = defaultdict(list)
for ns, nl in product(n_samples, n_labels):
    local_result = []
    for i in range(repeat):
        print(f"Repetition {i+1} for n_samples={ns} and n_labels={nl}")
        x = rng.integers(low=0, high=nl, size=ns)
        y = rng.integers(low=0, high=nl, size=ns)

        start = time()
        adjusted_mutual_info_score(x, y)
        end = time()
        local_result.append(end - start)

    result["n_samples"].append(ns)
    result["n_labels"].append(nl)
    result["mean_time"].append(np.mean(local_result))

result = pd.DataFrame(result)
plt.figure("Adjusted Mutual Info Score Benchmarks against number of Labels")
for n_sample in n_samples:
    samples = result[result["n_samples"] == n_sample]
    plt.plot(
        samples["n_labels"], samples["mean_time"], label=f"{str(n_sample)} samples"
    )
plt.xlabel("n_labels")
plt.ylabel("Time (s)")
plt.legend()

plt.figure("Adjusted Mutual Info Score Benchmarks against number of Samples")
for n_label in n_labels:
    labels = result[result["n_labels"] == n_label]
    plt.plot(labels["n_samples"], labels["mean_time"], label=f"{str(n_label)} labels")
plt.xlabel("n_samples")
plt.ylabel("Time (s)")
plt.legend()

plt.show()
