from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_predict
from sklearn.datasets import make_classification
from sklearn.utils.mocking import CheckingClassifier, MockDataFrame

list_check = lambda x: isinstance(x, list)
clf = CheckingClassifier(check_X=list_check)

X, y = make_classification(n_classes=2, n_samples=10)
predictions = cross_val_predict(LogisticRegression(), X.tolist(), y.tolist(), method='decision_function')
