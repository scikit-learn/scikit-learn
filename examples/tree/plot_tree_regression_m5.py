"""
===================================================================
Decision Tree Regression
===================================================================

A 1D regression with M5P decision tree.

The :ref:`decision trees <tree>` is
used to fit a sine curve with addition noisy observation. As a result, it
learns local linear regressions approximating the sine curve.

We can see that if the maximum depth of the tree (controlled by the
`max_depth` parameter) is set too high, the decision trees learn too fine
details of the training data and learn from the noise, i.e. they overfit.
"""
print(__doc__)

# Import the necessary modules and libraries
import numpy as np
from sklearn.tree import M5Prime, export_text_m5
import matplotlib.pyplot as plt

# Create a random dataset
rng = np.random.RandomState(1)
X = np.sort(5 * rng.rand(80, 1), axis=0)
y = np.sin(X).ravel()
y[::5] += 0.5 * (0.5 - rng.rand(16))

# Fit regression model
regr_1 = M5Prime(use_smoothing=False, use_pruning=False)
regr_1_label = "Tree 1 (use_smoothing=False, use_pruning=False)"
regr_1.fit(X, y)
regr_2 = M5Prime(use_smoothing=False)
regr_2_label = "Tree 2 (use_smoothing=False, use_pruning=True)"
regr_2.fit(X, y)
regr_3 = M5Prime(smoothing_constant_ratio=0.001)
regr_3_label = "Tree 3 (smoothing_constant_ratio=0.001, use_pruning=True)"
regr_3.fit(X, y)

# Predict
X_test = np.arange(0.0, 5.0, 0.01)[:, np.newaxis]
y_1 = regr_1.predict(X_test)
y_2 = regr_2.predict(X_test)
y_3 = regr_3.predict(X_test)

# Print the trees
print("\n----- %s" % regr_1_label)
print(regr_1.as_pretty_text())
print("\n----- %s" % regr_2_label)
print(regr_2.as_pretty_text())
print("\n----- %s" % regr_3_label)
print(export_text_m5(regr_3, out_file=None))  # equivalent to as_pretty_text

# Plot the results
plt.figure()
plt.scatter(X, y, s=20, edgecolor="black",
            c="darkorange", label="data")
plt.plot(X_test, y_1, color="cornflowerblue", label=regr_1_label, linewidth=2)
plt.plot(X_test, y_2, color="yellowgreen", label=regr_2_label, linewidth=2)
plt.plot(X_test, y_3, color="green", label=regr_3_label, linewidth=2)
plt.xlabel("data")
plt.ylabel("target")
plt.title("Decision Tree Regression")
plt.legend()
plt.show()
