"""
====================
Group K-Fold methods
====================

This example demonstrates when the stratify option of GroupKFold has an
advantage.
"""
from matplotlib import pyplot as plt
import numpy as np
from sklearn.model_selection import GroupKFold

print(__doc__)

rng = np.random.RandomState(0)
n_samples = 1000
n_groups = 100
X = np.arange(n_samples)
y = np.sort(rng.normal(size=n_samples))
groups = np.sort(rng.randint(0, n_groups, n_samples))

fig, axes = plt.subplots(1, 3, figsize=(18, 4), sharex=True, sharey=True)
for n, method in enumerate(('balance', 'stratify', 'shuffle')):
    cv = GroupKFold(2, method=method)
    for m, (train, test) in enumerate(cv.split(X, y, groups)):
        axes[n].hist(y[test], bins=20, histtype='step')
        print('%s fold %d: %d items' % (method, m + 1, len(test)))
    axes[n].set_xlabel(method)

plt.show()
