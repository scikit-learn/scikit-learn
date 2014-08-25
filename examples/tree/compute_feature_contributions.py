"""
================================================================
Compute the exact feature contributions for each prediction 
================================================================

Compute exactly how much each feature contributes to the final 
prediction in the decision tree regressor.
This enables one to understand why the model is giving a certain
result.

"""

import sys 
from sklearn.tree import DecisionTreeRegressor 
from sklearn.datasets import load_boston
from collections import defaultdict
import numpy as np

data = load_boston()
dt = DecisionTreeRegressor()

#Compute fature contributins on three samples. Use the rest of the data
#as training set
train = data.data[:-3]
target = data.target[:-3]
test = data.data[-3:]

dt.fit(train, target)

prediction = dt.predict(test)
#compute the decision paths from root to leaf for each sample
paths = dt.decision_paths(test)

contribution_list = []
for path in paths:
    contributions = defaultdict(int)
    for i in range(len(path)):
        if i == len(path) - 1 or path[i+1] == -1:
            break
        node_id = path[i]
        next_node_id = path[i+1]
        #feature contribution at a node is the difference of the mean of the 
        #decision node, and the mean at the following child node
        contributions[dt.tree_.feature[node_id]] += \
            dt.tree_.value[next_node_id][0][0] - dt.tree_.value[node_id][0][0]
    contribution_list.append(contributions)
  
print "Training set mean:", np.mean(target)
print
for i, contribs in enumerate(contribution_list):
    print "Prediction %s: %s " %(i + 1, prediction[i])
    print "Top five contributing features:"
    clist = sorted(contribs.items(), key = lambda x: -abs(x[1]))[:5]
    for key, val in clist:
        print data.feature_names[key], val
    print
        
