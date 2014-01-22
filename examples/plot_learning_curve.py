"""
========================
Plotting Learning Curves
========================

A learning curve shows the validation and training score of a learning
algorithm for varying numbers of training samples. It is a tool to
find out how much we benefit from adding more training data. If both
the validation score and the training score converge to a value that is
too low, we will not benefit much from more training data and we will
probably have to use a learning algorithm or a parametrization of the
current learning algorithm that can learn more complex concepts (i.e.
has a lower bias).

In this example, on the left side the learning curve of a naive Bayes
classifier is shown for the digits dataset. Note that the training score
and the cross-validation score are both not very good at the end. However,
the shape of the curve can be found in more complex datasets very often:
the training score is very high at the beginning and decreases and the
cross-validation score is very low at the beginning and increases. On the
right side we see the learning curve of an SVM with RBF kernel. We can
see clearly that the training score is still around the maximum and the
validation score could be increased with more training samples.
"""
print(__doc__)

import matplotlib.pyplot as plt
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.datasets import load_digits
from sklearn.learning_curve import learning_curve


digits = load_digits()
X, y = digits.data, digits.target

plt.figure()
plt.title("Learning Curve (Naive Bayes)")
plt.xlabel("Training examples")
plt.ylabel("Score")
train_sizes, train_scores, test_scores = learning_curve(
    GaussianNB(), X, y, cv=10, n_jobs=1)
plt.plot(train_sizes, train_scores, label="Training score")
plt.plot(train_sizes, test_scores, label="Cross-validation score")
plt.legend(loc="best")

plt.figure()
plt.title("Learning Curve (SVM, RBF kernel, $\gamma=0.001$)")
plt.xlabel("Training examples")
plt.ylabel("Score")
train_sizes, train_scores, test_scores = learning_curve(
    SVC(gamma=0.001), X, y, cv=10, n_jobs=1)
plt.plot(train_sizes, train_scores, label="Training score")
plt.plot(train_sizes, test_scores, label="Cross-validation score")
plt.legend(loc="best")

plt.show()
