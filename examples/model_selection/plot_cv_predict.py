"""
====================================
Plotting Cross-Validated Predictions
====================================

This example shows how to use
:func:`~sklearn.model_selection.cross_val_predict` to visualize prediction
errors.
"""
import matplotlib.pyplot as plt
from sklearn.datasets import load_diabetes
from sklearn.linear_model import LinearRegression
from sklearn.metrics import PredictionErrorDisplay
from sklearn.model_selection import cross_val_predict

X, y = load_diabetes(return_X_y=True)
lr = LinearRegression()

# cross_val_predict returns an array of the same size as `y` where each entry
# is a prediction obtained by cross validation:
y_pred = cross_val_predict(lr, X, y, cv=10)

display = PredictionErrorDisplay(y_true=y, y_pred=y_pred)
display.plot()
plt.show()
