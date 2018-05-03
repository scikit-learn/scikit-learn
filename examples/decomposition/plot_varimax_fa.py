"""
===============================================================
Model selection with Probabilistic PCA and Factor Analysis (FA)
===============================================================

Probabilistic PCA and Factor Analysis are probabilistic models.
The consequence is that the likelihood of new data can be used
for model selection and covariance estimation.
Here we compare PCA and FA with cross-validation on low rank data corrupted
with homoscedastic noise (noise variance
is the same for each feature) or heteroscedastic noise (noise variance
is the different for each feature). In a second step we compare the model
likelihood to the likelihoods obtained from shrinkage covariance estimators.

One can observe that with homoscedastic noise both FA and PCA succeed
in recovering the size of the low rank subspace. The likelihood with PCA
is higher than FA in this case. However PCA fails and overestimates
the rank when heteroscedastic noise is present. Under appropriate
circumstances the low rank models are more likely than shrinkage models.

The automatic estimation from
Automatic Choice of Dimensionality for PCA. NIPS 2000: 598-604
by Thomas P. Minka is also compared.

"""

# Authors: Jona Sassenhagen
# License: BSD 3 clause

import matplotlib.pyplot as plt
import numpy as np

from sklearn.decomposition import FactorAnalysis
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import load_iris

print(__doc__)

# #############################################################################
# Load Iris data
data = load_iris()
X = StandardScaler().fit_transform(data["data"])
feature_names = data["feature_names"]

n_comps = 2

# #############################################################################
# Run factor analysis with Varimax rotation
fig, axes = plt.subplots(ncols=2, figsize=(10, 8))

for ax, method in zip(axes, (None, "varimax")):
    fa = FactorAnalysis(n_components=n_comps, rotation=method).fit(X)

    components = fa.components_.T
    print(method, ":\n", components, end="\n\n")

    vmax = np.abs(components).max()
    ax.imshow(components, cmap="RdBu_r", vmax=vmax, vmin=-vmax)
    ax.set_yticks(np.arange(len(feature_names)))
    if ax.is_first_col():
        ax.set_yticklabels(feature_names)
    else:
        ax.set_yticklabels([])
    ax.set_xlabel("Rotation: " + str(method))
fig.suptitle("Factors")
plt.show()
