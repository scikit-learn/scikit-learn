r"""
========================================================
Post pruning decision trees with cost complexity pruning
========================================================

In this example, decision tree classifiers are trained with a post pruning
technique called minimal cost complexity pruning. This technique is
parameterized by the complexity parameter, ``ccp_alpha``. Greater values of
``ccp_alpha`` will prune more of the tree, thus creating smaller trees.
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
from sklearn.tree import DecisionTreeClassifier, plot_tree

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
plt.show()

###############################################################################
# Overview of cost complexity pruning
# -----------------------------------
# In this section, we will go through one iteration of the cost complexity
# pruning algorithm.
#
# Train a small decision tree with a maximum depth of 3
# .....................................................
clf_orig = DecisionTreeClassifier(random_state=0, max_depth=3)
clf_orig.fit(X_train, y_train)

###############################################################################
# Calcuate cost complexity measure
# ............................................
# Bubble up the leaf's :math:`R(t)` to its ancestors to calcuate,
# :math:`R(T_t)`.
weighted_n_samples = clf_orig.tree_.weighted_n_node_samples
r_t = clf_orig.tree_.impurity * weighted_n_samples / weighted_n_samples[0]
node_count = clf_orig.tree_.node_count

# Get node_id of parents
parents = np.zeros(node_count, dtype=int)
# leaf boolean mask
is_leaf = np.zeros(node_count, dtype=bool)
# number of leaves
n_leaves = np.zeros(node_count, dtype=int)

child_l = clf_orig.tree_.children_left
child_r = clf_orig.tree_.children_right

stack = [(0, -1)]
while stack:
    node_id, parent = stack.pop()
    parents[node_id] = parent
    if child_l[node_id] != -1:
        stack.append((child_l[node_id], node_id))
        stack.append((child_r[node_id], node_id))
    else:
        is_leaf[node_id] = True

# calcuate cost complexity measure
r_T = np.zeros(node_count, dtype=float)
leaf_indices = np.where(is_leaf)[0]
for leaf in leaf_indices:
    # bubble up leave values
    stack = [leaf]
    leaf_impurity = r_t[leaf]
    while stack:
        node_id = stack.pop()
        r_T[node_id] += leaf_impurity
        if parents[node_id] != -1:
            n_leaves[parents[node_id]] += 1
            stack.append(parents[node_id])

###############################################################################
#
# Find weakest link
# .................
# The weakest link is the node that has the smallest :math:`\alpha_{eff}`
alpha = (r_t - r_T) / (n_leaves - 1)

# leaves can not be pruned
alpha[is_leaf] = np.inf

weakest_node_id = np.argmin(alpha)
min_alpha = alpha[weakest_node_id]
n_samples = clf_orig.tree_.n_node_samples[weakest_node_id]
print("weakest node id: {}, alpha: {}, samples: {}".format(
      weakest_node_id, min_alpha, n_samples))


###############################################################################
#
# Plotting tree with node weakest link
# ....................................
# The weakest link has a background color in sky blue.
ax = plot_tree(clf_orig)
ax[weakest_node_id].set_backgroundcolor("skyblue")
plt.show()

###############################################################################
#
# Train a decision tree with pruning
# ..................................
# The cost complexity pruning algorithm prunes nodes till the minimal
# effective alpha is greater than ``ccp_alpha``. The resulting tree
# has the leaves of the weakest link pruned.
clf_pruned = DecisionTreeClassifier(random_state=0, max_depth=3,
                                    ccp_alpha=min_alpha + 0.0001)
clf_pruned.fit(X_train, y_train)
ax = plot_tree(clf_pruned)
ax[weakest_node_id].set_backgroundcolor("skyblue")
plt.show()
