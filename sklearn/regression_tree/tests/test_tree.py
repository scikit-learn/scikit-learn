from sklearn.datasets import load_diabetes
from sklearn.datasets import load_boston
from sklearn.regression_tree import RegressionTree
from sklearn.tree.tree import DecisionTreeRegressor
from sklearn.tree import export_graphviz

# data = load_diabetes()
data = load_boston()
X, y = data.data, data.target
X, y = X[:20], y[:20]

reg_tree = RegressionTree(random_state=42)
reg_tree.fit(X, y)
reg_tree.predict(X)
export_graphviz(reg_tree, out_file='reg_tree.dot')
print('Score regression tree: {}'.format(reg_tree.score(X, y)))

dec_tree = DecisionTreeRegressor(random_state=42)
dec_tree.fit(X, y)
dec_tree.predict(X)
export_graphviz(dec_tree, out_file='dec_tree.dot')
print('Score regression tree: {}'.format(dec_tree.score(X, y)))
