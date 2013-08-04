from sklearn.datasets import load_digits
from sklearn.neural_network import MultilayerPerceptronClassifier
from itertools import cycle, izip

digits = load_digits(n_class=2)
X, y = digits.data, digits.target
print 'Digits Dataset loaded'
for algorithm in ['sgd', 'l-bfgs']:
    for activation in ['tanh', 'logistic']:
        clf = MultilayerPerceptronClassifier(
            algorithm=algorithm,
            activation=activation,
            random_state=1).fit(X, y)
        print(
            "training accuracy for %s-based %s: %f" %
            (activation, algorithm, clf.score(X, y)))
