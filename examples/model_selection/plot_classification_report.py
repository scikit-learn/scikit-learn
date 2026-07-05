"""
=================================================
Visualizing a classification report as a heatmap
=================================================

The :func:`~sklearn.metrics.classification_report` function produces a text
summary of the main classification metrics (precision, recall and f1-score) for
each class. When there are more than a handful of classes, or when comparing
several models, the text form can be hard to scan.

:class:`~sklearn.metrics.ClassificationReportDisplay` renders the same numbers
as a heatmap. The precision, recall and f1-score columns share a fixed ``0``-``1``
color scale, so the color of a cell reflects *classification quality* directly --
unlike a raw confusion matrix, where color tracks the number of samples. The
``support`` column (the number of true instances of each class) is shown as
plain text since it is a count, not a score.

This example uses the digits dataset and a support vector classifier.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import matplotlib.pyplot as plt

from sklearn.datasets import load_digits
from sklearn.metrics import ClassificationReportDisplay, classification_report
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC

# %%
# Train a classifier
# ------------------
#
# We split the digits dataset and fit a support vector classifier.
X, y = load_digits(return_X_y=True)
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

# A slightly under-tuned ``gamma`` leaves a few classes imperfect, so the
# heatmap has some color variation to look at.
classifier = SVC(gamma=0.005).fit(X_train, y_train)

# %%
# Plot from an estimator
# ----------------------
#
# The quickest way to obtain the display is
# :meth:`~sklearn.metrics.ClassificationReportDisplay.from_estimator`, which
# computes the report from the fitted estimator and the test data.
disp = ClassificationReportDisplay.from_estimator(classifier, X_test, y_test)
disp.ax_.set_title("Classification report on the digits test set")
plt.show()

# %%
# Reading the plot
# ----------------
#
# Each row is a class (the ten digits), followed by the ``macro avg`` and
# ``weighted avg`` summary rows. Cells closer to ``1.0`` are brighter (better);
# darker cells highlight classes where the model is weaker -- here the precision
# of one digit lags behind the rest, which is immediately visible as a darker
# cell. Because the color scale is fixed to ``[0, 1]``, a class with few samples
# but a poor score stands out just as clearly as a well-represented one.

# %%
# Plot from an existing report
# ----------------------------
#
# If you already have a report dictionary from
# :func:`~sklearn.metrics.classification_report`, pass it straight to the
# constructor. This is handy when the report was computed elsewhere, or when you
# want to render several stored reports side by side.
report = classification_report(y_test, classifier.predict(X_test), output_dict=True)

ClassificationReportDisplay(report).plot()
plt.show()
