import sys 
sys.path.insert(0, 'c:/projects/scikit-learn/build/lib.win32-2.7/')
from sklearn.tree import DecisionTreeRegressor 
from sklearn.datasets import load_boston
from collections import defaultdict
import numpy as np

data = load_boston()
dt = DecisionTreeRegressor()

train = data.data[:-3]
target = data.target[:-3]

test = data.data[-3:]

dt.fit(train, target)
prediction = dt.predict(test)

paths = dt.predict(test, return_paths = True)
contribution_list = []
for path in paths:
    base = 0
    pred = 0
    contributions = defaultdict(int)
    for i in range(len(path)):
        if path[i+1] == -1:
            break
        node_id = path[i]
        next_node_id = path[i+1]
        
        contributions[dt.tree_.feature[node_id]] += dt.tree_.value[next_node_id][0][0] - dt.tree_.value[node_id][0][0]
        
        
    contribution_list.append(contributions)
  
print "Training set mean:", np.mean(target)
print
for i, contribs in enumerate(contribution_list):
    print "Prediction %s: %s " %(i + 1, prediction[i])
    print "Top three contributing features:"
    clist = sorted(contribs.items(), key = lambda x: -abs(x[1]))[:3]
    for key, val in clist:
        print data.feature_names[key], val
    print
        
