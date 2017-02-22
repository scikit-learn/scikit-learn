from sklearn.datasets import load_diabetes
from sklearn.regression_tree import RegressionTree
from sklearn.tree.tree import DecisionTreeRegressor
from sklearn.tree import export_graphviz

data = load_diabetes()
X, y = data.data, data.target
X, y = X[:20], y[:20]

reg_tree = RegressionTree()
reg_tree.fit(X, y)
reg_tree.predict(X)
export_graphviz(reg_tree, out_file='reg_tree.dot')
