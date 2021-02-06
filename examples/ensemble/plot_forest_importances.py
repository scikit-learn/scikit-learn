"""
=========================================
Feature importances with forests of trees
=========================================

This examples shows the use of forests of trees to evaluate the importance of
features on an artificial classification task. The red bars are
the impurity-based feature importances of the forest,
along with their inter-trees variability.

As expected, the plot suggests that 3 features are informative, while the
remaining are not.

.. warning::
    Impurity-based feature importances can be misleading for high cardinality
    features (many unique values). 
"""
print(__doc__)

# %%
import time
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from sklearn.datasets import make_classification

# Build a classification task using 3 informative features
X, y = make_classification(n_samples=1000,
                           n_features=10,
                           n_informative=3,
                           n_redundant=0,
                           n_repeated=0,
                           n_classes=2,
                           random_state=0,
                           shuffle=False)


# %%
# Feature permutation importance on full model
# --------------------------------------------
# Another way to calculate the feature importances is by 
# using permutation_importance from the inspection submodule.
# Here, we'll fit an Extra-Trees Classifier to the data and 
# calculate the feature importances on the full data set.
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.inspection import permutation_importance

feature_names = [f'feature {i}' for i in range(X.shape[1])]

forest = ExtraTreesClassifier(n_estimators=250,
                              random_state=42)

start_time = time.time()                        
forest.fit(X, y)
result = permutation_importance(forest, X, y, n_repeats=10,
                                random_state=42, n_jobs=2)
forest_importances = pd.Series(
    result.importances_mean, index=feature_names)
elapsed_time = time.time() - start_time
print(f"Elapsed time to fit and compute the importances: "
      f"{elapsed_time:.3f} seconds")

# The computation for full permutation importance is more costly.
# Features are shuffled n times and the model refitted to estimate
# the importance of it. Please see :ref:`permutation_importance` for
# more details.

# We can now plot the importance ranking.

ax = forest_importances.plot.bar(yerr=result.importances_std)
ax.set_title("Feature importances using permutation on full model")
ax.set_ylabel("Mean accuracy decrease")
plt.show()

# Same as in the OOB feature importance three most important features
# are detected.
