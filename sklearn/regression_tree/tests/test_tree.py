from sklearn.datasets import load_diabetes
from sklearn.regression_tree import RegressionTree
from sklearn.tree.tree import DecisionTreeRegressor
from sklearn.tree import export_graphviz

data = load_diabetes()

reg_tree = RegressionTree(max_depth=4)
reg_tree.fit(data.data, data.target)
reg_tree.predict(data.data)
export_graphviz(reg_tree, out_file='reg_tree.dot')

dec_tree = DecisionTreeRegressor(max_depth=4)
dec_tree.fit(data.data, data.target)
export_graphviz(dec_tree, out_file='dec_tree.dot')
