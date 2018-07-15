from sklearn.linear_model import LogisticRegression
from sklearn import datasets
from sklearn.model_selection import GridSearchCV


digits = datasets.load_digits()
iris = datasets.load_iris()
X = iris.data[:, [2, 3]]
y = iris.target


model3 = LogisticRegression(verbose=101, solver='lbfgs', max_iter=101)
model3.fit(X, y)


if __name__ == '__main__':
    print('Hello World')