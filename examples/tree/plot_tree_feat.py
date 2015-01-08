"""
===================================================================
Decision Tree Feature Extraction
===================================================================

Obtaining features from decision trees.

A dataset can be transformed using a decision tree's apply() method
in two ways:

1) Reducing the number of classes to predict. By selecting max_leaf_nodes
to be a value less than the total number of classes in a classification
problem, one can obtain a dataset with a reduced number of classes.

2) Creating a new sparse feature representation of the data.
Each sample will be transformed to a vector of the size of the number of 
leafs in the decision tree. Each leaf is assigned an index in this vector.
If the sample falls into a given leaf, the value at that leaf's index in the
vector is 1; otherwise it is 0. Only one value in the array can be 1.
This sparse, high-dimensional representation may be useful for increasing data
separability.

Note that in the below double bar graph demonstrating the first
transformation, all 
setosas fall into the first leaf, and none into the second leaf. Similarly, 
all versicolors and virginicas fall only into the second leaf. This suggests
that virginicas and versicolors are more similar to each other than to 
setosas.

"""
print(__doc__)

import numpy as np
import matplotlib.pyplot as plt

from sklearn import tree
from sklearn import ensemble
from sklearn.datasets import load_iris

max_leaves = 2
iris = load_iris()
clf = tree.DecisionTreeClassifier(max_leaf_nodes=max_leaves)
X = iris['data']
y = iris['target']

clf.fit(X,y)
#1
y_reduced = clf.apply(X) #Now only two classes instead of three.

bar_width = .35
opacity = 0.4
index = np.arange(3)

leaf_class_colors = {}
leaf_class_colors.update(zip(range(np.max(y_reduced)), ['r', 'b']))

new_classes = []
for i in (1,2):
    new_classes.append(np.array([np.sum(y[y_reduced == i] == 0), \
                                        np.sum(y[y_reduced == i] == 1), \
                                        np.sum(y[y_reduced == i] == 2) \
                                    ]))

for i in xrange(np.max(y_reduced)):
    plt.bar(index + i * bar_width, new_classes[i], bar_width, alpha=opacity, \
            color=leaf_class_colors[i], label="Leaf " +str(i + 1))

plt.title("The assignment of each original class to new leaf index classes")
plt.xticks(index + bar_width, iris['target_names'])
plt.xlabel("Original class")
plt.ylabel("Number in each new leaf class")
plt.legend()


#2
# We don't need to use a decision tree with a constrained number of leaves,
# but we do so here for the convenience of using the same classifier to
# demonstrate part 1.
X_trans = np.zeros((y_reduced.size, max_leaves))
for i in xrange(max_leaves):
    X_trans[:,i] = y_reduced == i + 1 #Add 1 because leaf indexing begins at 1

#For the graph to be built in the documentation, plt.show() must be called last.
plt.show()
