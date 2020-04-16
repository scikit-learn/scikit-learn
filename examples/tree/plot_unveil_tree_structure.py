"""
=========================================
Understanding the decision tree structure
=========================================

The decision tree structure can be analysed to gain further insight on the
relation between the features and the target to predict. In this example, we
show how to retrieve:

- the binary tree structure;
- the depth of each node and whether or not it's a leaf;
- the nodes that were reached by a sample using the ``decision_path`` method;
- the leaf that was reached by a sample using the apply method;
- the rules that were used to predict a sample;
- the decision path shared by a group of samples.

"""
import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.datasets import load_iris
from sklearn.tree import DecisionTreeClassifier
from sklearn import tree

##############################################################################
# Train tree classifier
# ---------------------
# First we fit a :class:`tree.DecisionTreeClassifier` using the
# :func:`datasets.load_iris` dataset.

iris = load_iris()
X = iris.data
y = iris.target
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

clf = DecisionTreeClassifier(max_leaf_nodes=3, random_state=0)
clf.fit(X_train, y_train)

##############################################################################
# Tree structure
# --------------
#
# The decision clf has an attribute called ``tree_`` which allows access
# to low level attributes such as ``node_count``, the total number of nodes,
# and ``max_depth``, the maximal depth of the tree. It also stores the
# entire binary tree structure, represented as a number of parallel arrays. The
# i-th element of each array holds information about the node i. Node 0 is the
# tree's root. Some of the arrays only apply to either leaves or split nodes.
# In this case the values of the nodes of the other type is arbitrary. For
# example, the attributes ``feature`` and ``threshold`` only apply to split
# nodes. The values corresponding to leaf nodes in these arrays are therefore
# arbitrary.
#
# Among those arrays, we have:
#   - ``children_left`` - id of the left child of node i, -1 if leaf node
#   - ``children_right`` - id of the right child of node i, -1 if leaf node
#   - ``feature`` - feature used for splitting node i
#   - ``threshold`` - threshold value at node i
#   - ``n_node_samples`` - the number of of training samples reaching node i
#   - ``impurity`` - the impurity at node i
#
# Using the arrays, we can traverse the tree structure to compute various
# properties. Below, we will compute the depth of each node and whether or not
# it is a leaf.

n_nodes = clf.tree_.node_count
children_left = clf.tree_.children_left
children_right = clf.tree_.children_right
feature = clf.tree_.feature
threshold = clf.tree_.threshold

node_depth = np.zeros(shape=n_nodes, dtype=np.int64)
is_leaves = np.zeros(shape=n_nodes, dtype=bool)
stack = [(0, -1)]  # start with the root node id (0) and its parent depth (-1)
while len(stack) > 0:
    # `pop` ensures each node is only visited once
    node_id, parent_depth = stack.pop()
    # Increment parent depth to get current node depth
    node_depth[node_id] = parent_depth + 1

    # If the left and right child of a node is not the same we have a split
    # node
    is_split_node = children_left[node_id] != children_right[node_id]
    # If a split node, append left and right children and depth to `stack`
    # so we can loop through them
    if is_split_node:
        stack.append((children_left[node_id], parent_depth + 1))
        stack.append((children_right[node_id], parent_depth + 1))
    else:
        is_leaves[node_id] = True

print("The binary tree structure has %s nodes and has "
      "the following tree structure:\n"
      % n_nodes)
for i in range(n_nodes):
    if is_leaves[i]:
        print("%snode=%s is a leaf node." % (node_depth[i] * "\t", i))
    else:
        print("%snode=%s is a split node: "
              "go to node %s if X[:, %s] <= %s else to node %s."
              % (node_depth[i] * "\t",
                 i,
                 children_left[i],
                 feature[i],
                 threshold[i],
                 children_right[i],
                 ))
print()

##############################################################################
# We can compare the above output to the plot of the decision tree.

tree.plot_tree(clf)

##############################################################################
# Decision path
# -------------
#
# We can also retrieve the decision path of each sample. The ``decision_path``
# method allows us to retrieve the nodes samples go through. A non zero element
# in the indicator matrix at position (i, j) indicates that the sample i goes
# through the node j. For one sample, i, the positions of the non zero
# element(s) in row i of the indicator matrix indicate the node(s) that sample
# goes through.
#
# The leaf ids reached by each sample can be obtained with the ``apply``
# method. We can use this and the ``decision_path`` to  get the tests that were
# used to predict a sample or a group of samples. First, let's do it for one
# sample. Note: ``node_index`` is a sparse matrix.

node_indicator = clf.decision_path(X_test)
leaf_id = clf.apply(X_test)

sample_id = 0
# obtain ids of the nodes `sample_id` goes through
node_index = node_indicator.indices[node_indicator.indptr[sample_id]:
                                    node_indicator.indptr[sample_id + 1]]

print('Rules used to predict sample %s:\n' % sample_id)
for node_id in node_index:
    # continue to the next node if it is a leaf node
    if leaf_id[sample_id] == node_id:
        continue

    if (X_test[sample_id, feature[node_id]] <= threshold[node_id]):
        threshold_sign = "<="
    else:
        threshold_sign = ">"

    print("decision id node %s : (X_test[%s, %s] = %s) %s %s)"
          % (node_id,
             sample_id,
             feature[node_id],
             X_test[sample_id, feature[node_id]],
             threshold_sign,
             threshold[node_id]))

##############################################################################
# For a group of samples, we can determine common nodes the samples go through.

sample_ids = [0, 1]
# boolean array indicating the nodes both samples go through
common_nodes = (node_indicator.toarray()[sample_ids].sum(axis=0) ==
                len(sample_ids))
# position in array equals node id
common_node_id = np.arange(n_nodes)[common_nodes]

print("\nThe following samples %s share the node(s) %s in the tree"
      % (sample_ids, common_node_id))
print("This is %s %% of all nodes." % (100 * len(common_node_id) / n_nodes,))
