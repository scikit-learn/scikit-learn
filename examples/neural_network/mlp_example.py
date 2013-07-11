from sklearn.datasets import load_digits
from sklearn.neural_network import MLPClassifier
from itertools import cycle, izip

digits = load_digits()
X, y = digits.data, digits.target
print 'Digits Dataset loaded'
for algorithm in ['sgd', 'l-bfgs']:
    for activation in ['tanh', 'logistic']:
        clf = MLPClassifier(
            algorithm=algorithm,
            activation=activation,
            random_state=1).fit(X, y)
        print(
            "training accuracy for %s-based %s: %f" %
            (activation, algorithm, clf.score(X, y)))
