"""
========================================================
Post pruning decision trees with cost complexity pruning
========================================================

In this example, decision tree classifiers are trained with a post pruning
technique called minimal cost complexity pruning. This technique is
parameterized by the complexity parameter, ``ccp_alpha``. Greater values of
``ccp_alpha`` will prune more of the tree, thus creating smaller trees.
"""

print(__doc__)
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.datasets import load_breast_cancer
from sklearn.tree import DecisionTreeClassifier

###############################################################################
# Total impurity of leaves vs effective alphas of pruned tree
# ---------------------------------------------------------------
# Minimal cost complexity pruning recursively finds the node with the
# "weakest link". The weakest link is characterized by an effective alpha,
# where the nodes with the smallest effective alpha are pruned first.
# :func:`~sklearn.tree.DecisionTreeClassifier.cost_complexity_pruning_path`
# returns the effective alphas and the corresponding total leaf impurities
# at each step of the pruning process. With higher effective alphas, more
# of the tree is pruned, which increases the total impurity of its leaves.
X, y = load_breast_cancer(return_X_y=True)
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

clf = DecisionTreeClassifier(random_state=0)
path = clf.cost_complexity_pruning_path(X_train, y_train)
ccp_alphas, impurities = path.ccp_alphas, path.impurities
fig, ax = plt.subplots()
ax.plot(ccp_alphas, impurities, marker='o')
ax.set_xlabel("effective alpha")
ax.set_ylabel("total impurity of leaves")
ax.set_title("Total Impurity vs effective alpha for training set")
plt.show()

###############################################################################
# Next, we train a decision tree using the effective alphas. The last element
# in ``clfs`` corresponds to a tree pruned till it has one node.
clfs = []
for ccp_alpha in ccp_alphas:
    clf = DecisionTreeClassifier(random_state=0, ccp_alpha=ccp_alpha)
    clf.fit(X_train, y_train)
    clfs.append(clf)
print("Number of nodes in the last tree is: {}".format(
      clfs[-1].tree_.node_count))

###############################################################################
# As alpha increases, more of the tree is pruned, resulting in trees with
# fewer nodes and less depth:
node_counts = [clf.tree_.node_count for clf in clfs]
fig, ax = plt.subplots()
ax.set_xlabel("alpha")
ax.set_ylabel("number of nodes")
ax.set_title("Number of nodes vs alpha")
ax.plot(ccp_alphas, node_counts, marker="o")
plt.show()

###############################################################################
# Accuracy vs alpha for training and testing sets
# ----------------------------------------------------
# We plot the training and testing accuracy as more of the tree is pruned. The
# last element in ``clfs`` and ``ccp_alphas`` is removed, because it is a tree
# with only one node.
clfs = clfs[:-1]
ccp_alphas = ccp_alphas[:-1]
train_scores = [clf.score(X_train, y_train) for clf in clfs]
test_scores = [clf.score(X_test, y_test) for clf in clfs]

fig, ax = plt.subplots()
ax.set_xlabel("alpha")
ax.set_ylabel("accuracy")
ax.set_title("Accuracy vs alpha for training and testing sets")
ax.plot(ccp_alphas, train_scores, label="train", marker="o")
ax.plot(ccp_alphas, test_scores, label="test", marker="o")
ax.legend()
plt.show()
