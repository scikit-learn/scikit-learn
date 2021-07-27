from sklearn.datasets import make_classification
X, y = make_classification()

from sklearn.linear_model import LogisticRegression
lr = LogisticRegression(C='wrong-value')

from sklearn.model_selection import cross_val_score
cross_val_score(lr, X, y, n_jobs=2)
