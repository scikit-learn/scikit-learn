"""
========================================================
Post pruning decision trees with cost complexity pruning
========================================================

In this example, decision tree classifiers are trained with a post pruning
technique called minimal cost complexity pruning. This technique is
parameterized by the complexity parameter, :math:`\alpha`. Greater values of
:math:`\alpha` will prune more of the tree, thus creating smaller trees.
"""
print(__doc__)

###############################################################################
# Train decision tree classifiers
# -------------------------------
# Train 40 decision tree classifiers with :math:`\alpha` from 0.00 to
# 0.40.
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.datasets import load_breast_cancer
from sklearn.tree import DecisionTreeClassifier

X, y = load_breast_cancer(return_X_y=True)
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

ccp_alphas = np.linspace(0, 0.04, 40)
clfs = []
for ccp_alpha in ccp_alphas:
    clf = DecisionTreeClassifier(random_state=0, ccp_alpha=ccp_alpha)
    clf.fit(X_train, y_train)
    clfs.append(clf)

###############################################################################
# Plot training and test scores vs alpha
# --------------------------------------
# Calculate and plot the training scores and test accuracy scores
# for our classifiers. With :math:`\alpha` equal to 0.0, the decision tree is
# overfitting with a 1.0 training accuracy score. As the decision tree is
# pruned the testing accuracy score increases up to a point and then decreases.
train_scores = []
test_scores = []
for clf in clfs:
    train_scores.append(clf.score(X_train, y_train))
    test_scores.append(clf.score(X_test, y_test))

fig, ax = plt.subplots()
ax.set_xlabel("alpha")
ax.set_ylabel("accuracy")
ax.set_title("Accuracy vs alpha for training and testing sets")
ax.plot(ccp_alphas, train_scores, label="train")
ax.plot(ccp_alphas, test_scores, label="test")
ax.legend()

###############################################################################
# Plot total number of nodes vs alpha
# -----------------------------------
# Plot the total number of nodes for our classifiers. As :math:`\alpha`
# increases, the number of nodes decreases.
node_counts = [clf.tree_.node_count for clf in clfs]
fig, ax = plt.subplots()
ax.set_xlabel("alpha")
ax.set_ylabel("number of nodes")
ax.set_title("Number of nodes vs alpha")
ax.plot(ccp_alphas, node_counts)
fig.show()
