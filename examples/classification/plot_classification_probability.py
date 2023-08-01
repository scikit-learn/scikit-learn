"""
===============================
Plot classification probability
===============================

In this example, we will plot the classification probability for different classifiers.

We will employ the following classifiers for classification:

1. Support Vector Classifier (SVC)
2. L1 penalized logistic regression
3. L2 penalized logistic regression

Each logistic regression will be evaluated in both the One-Vs-Rest and multinomial settings.
Additionally, we will employ Gaussian process classification.

We will work with the Iris dataset which consists four numerical features
(sepal length, sepal width, petal length, and petal width), where the
flowers class (setosa, versicolor or virginica) is the target.

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD 3 clause
"""

import matplotlib.pyplot as plt
import numpy as np

from sklearn import datasets
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.svm import SVC

# %%
# Load Iris dataset
# -------------------------
# First, we load the Iris data as a pandas dataframe.
iris = datasets.load_iris()
# We only take the first two features for visualization.
X = iris.data[:, 0:2]
y = iris.target

n_features = X.shape[1]

# %%
# CLASSIFIER INITIALIZATION
# -------------------------------------------------------------

# We will initialize five different classifiers:
#-L1 penalized logistic regression with multinomial setting.
#-L2 penalized logistic regression with multinomial setting.
#-L2 penalized logistic regression with One-Vs-Rest setting.
#This classifier is not inherently designed for multiclass classification. As a consequence,
#it faces more challenges in distinguishing between class 2 and class 3 compared to the other estimators.
#-Linear Support Vector Classifier (SVC) with C as the regularization parameter
#and probability set to True for obtaining classification probabilities.
#-Gaussian Process Classifier (GPC) using the specified Gaussian Radial Basis Function (RBF) kernel.


# Set the value of C for regularization in the logistic regression and Linear SVC classifiers.
C = 10

# Define the kernel for Gaussian Process Classifier (GPC).
kernel = 1.0 * RBF([1.0, 1.0])  # Gaussian Radial Basis Function (RBF) kernel for GPC

# Create different classifiers with their respective hyperparameters.
classifiers = {
    "L1 logistic": LogisticRegression(
        C=C, penalty="l1", solver="saga", multi_class="multinomial", max_iter=10000
    ),
    "L2 logistic (Multinomial)": LogisticRegression(
        C=C, penalty="l2", solver="saga", multi_class="multinomial", max_iter=10000
    ),
    "L2 logistic (OvR)": LogisticRegression(
        C=C, penalty="l2", solver="saga", multi_class="ovr", max_iter=10000
    ),
    "Linear SVC": SVC(kernel="linear", C=C, probability=True, random_state=0),
    "GPC": GaussianProcessClassifier(kernel),
}


# %%
# Plot classification probabilities for each classifier
# -----------------------------------------------------------
# In this section, we train multiple classifiers and then display their decision boundaries
# and classification probabilities to compare their performance on the Iris data.

n_classifiers = len(classifiers)

# Set the size and adjust the subplot positions of the figure
plt.figure(figsize=(3 * 2, n_classifiers * 2))
plt.subplots_adjust(bottom=0.2, top=0.95)

# Create a meshgrid for the decision boundary visualization
xx = np.linspace(3, 9, 100)
yy = np.linspace(1, 5, 100).T
xx, yy = np.meshgrid(xx, yy)
Xfull = np.c_[xx.ravel(), yy.ravel()]

# Loop through each classifier in the 'classifiers' dictionary and visualize the results
for index, (name, classifier) in enumerate(classifiers.items()):
    # Fit the classifier using the training data (X, y)
    classifier.fit(X, y)

    # Make predictions on the training data
    y_pred = classifier.predict(X)

    # Calculate the accuracy of the classifier on the training data
    accuracy = accuracy_score(y, y_pred)
    print("Accuracy (train) for %s: %0.1f%% " % (name, accuracy * 100))

    # View classification probabilities:
    probas = classifier.predict_proba(Xfull)
    n_classes = np.unique(y_pred).size

    # Create subplots for each class
    for k in range(n_classes):
        plt.subplot(n_classifiers, n_classes, index * n_classes + k + 1)
        plt.title("Class %d" % k)

        # Add ylabel to the leftmost subplots
        if k == 0:
            plt.ylabel(name)

        # Display the probabilities as an image
        imshow_handle = plt.imshow(
            probas[:, k].reshape((100, 100)), extent=(3, 9, 1, 5), origin="lower"
        )
        plt.xticks(())
        plt.yticks(())

        # Mark the actual data points on the plot if they belong to the current class
        idx = y_pred == k
        if idx.any():
            plt.scatter(X[idx, 0], X[idx, 1], marker="o", c="w", edgecolor="k")

# Create a color bar for the probability values
ax = plt.axes([0.15, 0.04, 0.7, 0.05])
plt.title("Probability")
plt.colorbar(imshow_handle, cax=ax, orientation="horizontal")

# Show the plot
plt.show()
