"""
=================================================
Pixel importances with a parallel forest of trees
=================================================

This example show the use of a forest of trees to evaluate
the importance of the pixels in an image classification task (digits)
based on impurity and permutation importance.
The hotter the pixel, the more important.

The code below also illustrates how the construction and the computation
of the predictions can be parallelized within multiple jobs.
"""
print(__doc__)

import matplotlib.pyplot as plt

# %%
# Loading the data and model fitting
# ----------------------------------
# We use the faces data from datasets submodules and split the dataset
# into training and testing subsets. Also, we'll set the number of cores
# to use for the tasks.
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split

# %%
# We select the number of cores to use to perform parallel fitting of
# the forest model. `-1` means use all available cores.
n_jobs = -1

# %%
# Load the faces dataset
data = load_digits()
X, y = data.data, data.target
X_train, X_test, y_train, y_test = train_test_split(
    X, y, stratify=y, random_state=42)

# %%
# A random forest classifier will be fitted to compute the feature importances.
from sklearn.ensemble import RandomForestClassifier

forest = RandomForestClassifier(
    n_estimators=750, n_jobs=n_jobs, random_state=42)

forest.fit(X_train, y_train)

# %%
# Feature importance based on mean decrease in impurity
# -----------------------------------------------------
# Feature importances are provided by the fitted attribute
# `feature_importances_` and they are computed as the mean and standard
# deviation of accumulation of the impurity decrease within each tree.
#
# .. warning::
#     Impurity-based feature importances can be misleading for high cardinality
#     features (many unique values). See :ref:`permutation_importance` as
#     an alternative below.
import time

start_time = time.time()
img_shape = data.images[0].shape
importances = forest.feature_importances_
elapsed_time = time.time() - start_time

print(f"Elapsed time to compute the importances: "
      f"{elapsed_time:.3f} seconds")

# %%
# Let's plot the impurity-based importance.
imp_reshaped = importances.reshape(img_shape)
plt.matshow(imp_reshaped, cmap=plt.cm.hot)
plt.title("Pixel importances using impurity values")
plt.colorbar()
plt.tight_layout()
plt.show()

# %%
# Feature importance based on feature permutation
# -----------------------------------------------
# Permutation feature importance overcomes limitations of the impurity-based
# feature importance: they do not have a bias toward high-cardinality features
# and can be computed on a left-out test set.
from sklearn.inspection import permutation_importance

start_time = time.time()
result = permutation_importance(
    forest, X_test, y_test, n_repeats=10,
    random_state=42, n_jobs=n_jobs)
elapsed_time = time.time() - start_time
print(f"Elapsed time to compute the importances: "
      f"{elapsed_time:.3f} seconds")

# %%
# The computation for full permutation importance is more costly. Features are
# shuffled n times and the model refitted to estimate the importance of it.
# Please see :ref:`permutation_importance` for more details. We can now plot
# the importance ranking.

plt.matshow(result.importances_mean.reshape(img_shape), cmap=plt.cm.hot)
plt.title("Pixel importances using permutation importance")
plt.colorbar()
plt.tight_layout()
plt.show()

# %%
# We can see similar areas are detected using both methods. Although
# the importances vary. We can see that permutation importance gives lower
# importance values on any single pixel, which matches the intuition:
# The class of a digit seen on an image depends on values of many pixels
# together rather than a few pixels.
