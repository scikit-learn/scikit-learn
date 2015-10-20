"""
=========================================
Understanding the decision tree structure
=========================================

The decision tree structure can be analysed to gain further insight on the
relation between the features and the target to predict. In this example, we
show how to retrieve:
    - the binary tree structure;
    - the nodes that were reached by a sample using the decision_paths method;
    - the leaf that was reached by a sample using the apply method;
    - the rules that were used to predict a sample;
    - the decision path shared by a group of samples.

"""
import numpy as np

from sklearn.cross_validation import train_test_split
from sklearn.datasets import load_iris
from sklearn.tree import DecisionTreeRegressor


iris = load_iris()
X = iris.data
y = iris.target
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

estimator = DecisionTreeRegressor(max_leaf_nodes=3, random_state=0)
estimator.fit(X_train, y_train)

# The decision estimator has an attribute called tree_  which stores the entire
# tree structure and allows access to low level attributes. The binary tree
# tree_ is represented as a number of parallel arrays. The i-th element of each
# array holds information about the node `i`. Node 0 is the tree's root. NOTE:
# Some of the arrays only apply to either leaves or split nodes, resp. In this
# case the values of nodes of the other type are arbitrary!
#
# Among those arrays, we have:
#   - left_child, id of the left child of the node
#   - right_child, id of the right child of the node
#   - feature, feature used for splitting the node
#   - threshold, threshold value at the node
#

# Using those arrays, we can parse the tree structure:

print("The binary tree structure has %s nodes and has "
      "the following tree structure:"
      % estimator.tree_.node_count)

for i in range(estimator.tree_.node_count):
    if estimator.tree_.children_left[i] == estimator.tree_.children_right[i]:
        print("node=%s leaf node." % i)
    else:
        print("node=%s test node: go to node %s if X[:, %s] <= %ss else to "
              "node %s."
              % (i,
                 estimator.tree_.children_left[i],
                 estimator.tree_.feature[i],
                 estimator.tree_.threshold[i],
                 estimator.tree_.children_right[i],
                 ))
print()

# First let's retrieve the decision path of each sample. The decision_paths
# method allows to retrieve the node indicator function. A non zero elements at
# position (i, j) indicates that the sample i goes sthrough the node j.

node_indicator = estimator.decision_paths(X_test)

# Similarly, we can also have the leaves ids reach by each sample.

leave_id = set(estimator.apply(X_test))

# Now, it's possible to get the tests that were used to predict a sample or
# a group of samples. First, let's make it for the  sample.

sample_id = 0
node_index = node_indicator.indices[node_indicator.indptr[sample_id]:
                                    node_indicator.indptr[sample_id + 1]]

print('Rules used to predict sample %s: ' % sample_id)
for i, node_id in enumerate(node_index):
    if node_id in leave_id:
        continue

    if (X_test[i, estimator.tree_.feature[node_id]] <=
            estimator.tree_.threshold[node_id]):
        threshold_sign = "<="
    else:
        threshold_sign = ">"

    print("rule %s from node %s : (X[%s, %s] (= %s) %s %s)"
          % (i,
             node_id,
             sample_id,
             estimator.tree_.feature[node_id],
             X_test[i, estimator.tree_.feature[node_id]],
             threshold_sign,
             estimator.tree_.threshold[node_id]))

# For a group of samples, we have the following common node.
sample_ids = [0, 1]
common_nodes = (node_indicator.toarray()[sample_ids].sum(axis=0) ==
                len(sample_ids))

common_node_id = np.arange(estimator.tree_.node_count)[common_nodes]

print("\nThe following sample %s shares the following path %s  in the tree"
      % (sample_ids, common_node_id))
print("It is %s %% of all nodes."
      % (len(common_node_id) / estimator.tree_.node_count * 100,))
