"""
=============================
Recursive feature elimination
=============================


This example demonstrates how :class:~sklearn.feature_selection.RFE can be used
to determine the importance of individual pixels when classifying handwritten digits.
RFE is a method that recursively removes the least significant features and retrains
the model, allowing us to rank features by their importance. the most important features
are assigned rank 1 and the higher ranking_, the less important. This ranking is also
encoded by the shades of blue. Ranking map is annotated with pixel numbers for clarity.
As expected, pixels in the center of the image are more predictive than pixels close to the edges.

.. note::

    See also :ref:`sphx_glr_auto_examples_feature_selection_plot_rfe_with_cross_validation.py`

"""  # noqa: E501

import matplotlib.pyplot as plt

from sklearn.datasets import load_digits
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

# Load the digits dataset
digits = load_digits()
X = digits.images.reshape((len(digits.images), -1))
y = digits.target

pipe = Pipeline(
    [
        ("scaler", StandardScaler()),
        ("rfe", RFE(estimator=LogisticRegression(), n_features_to_select=1, step=1)),
    ]
)

pipe.fit(X, y)
ranking = pipe.named_steps["rfe"].ranking_.reshape(digits.images[0].shape)

# Plot pixel ranking
plt.matshow(ranking, cmap=plt.cm.Blues)

# Add annotations for pixel numbers
for i in range(ranking.shape[0]):
    for j in range(ranking.shape[1]):
        plt.text(j, i, str(ranking[i, j]), ha="center", va="center", color="black")

plt.colorbar()
plt.title("Ranking of pixels with RFE\n(Logistic Regression)")
plt.show()
