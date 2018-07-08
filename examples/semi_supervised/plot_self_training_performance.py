import numpy as np
import matplotlib.pyplot as plt
from sklearn.semi_supervised.self_training import SelfTraining
from sklearn.utils import shuffle
from sklearn.svm import SVC
from sklearn.datasets import load_iris
from sklearn.model_selection import StratifiedKFold
from sklearn.semi_supervised.label_propagation import LabelPropagation
from sklearn.metrics import f1_score
from sklearn.base import clone

supervised_score = []
self_training_score = []
label_propagation_score = []
x_values = []

clf = SVC(probability=True, C=100, gamma=0.7, kernel='rbf')
self_training_clf = SelfTraining(
    clone(clf, safe=True), max_iter=100, threshold=0.8
)
ls = LabelPropagation()

for t in range(20, 80):
    x_values.append(t)
    X, y = load_iris(return_X_y=True)
    X, y = shuffle(X, y, random_state=42)
    y_true = y.copy()

    lim = t
    y[lim:] = -1

    supervised_score_temp = []
    self_training_score_temp = []
    label_propagation_score_temp = []

    skfolds = StratifiedKFold(n_splits=3, random_state=42)
    for train_index, test_index in skfolds.split(X, y):
        X_train = X[train_index]
        y_train = y[train_index]
        X_test = X[test_index]
        y_test = y[test_index]
        y_test_true = y_true[test_index]

        X_train_filtered = X_train[np.where(y_train != -1)]
        y_train_filtered = y_train[np.where(y_train != -1)]

        clf.fit(X_train_filtered, y_train_filtered)
        y_pred = clf.predict(X_test)
        supervised_score_temp.append(
            f1_score(y_test_true, y_pred, average='macro')
        )

        self_training_clf.fit(X_train, y_train)
        y_pred = self_training_clf.predict(X_test)
        self_training_score_temp.append(
            f1_score(y_test_true, y_pred, average='macro')
        )

        ls.fit(X_train, y_train)
        y_pred = ls.predict(X_test)
        label_propagation_score_temp.append(
            f1_score(y_test_true, y_pred, average='macro')
        )

    supervised_score.append(np.array(supervised_score_temp).mean())
    self_training_score.append(np.array(self_training_score_temp).mean())
    label_propagation_score.append(
        np.array(label_propagation_score_temp).mean()
    )


plt.figure(1)
plt.plot(x_values, supervised_score, label='Supervised')
plt.plot(x_values, self_training_score, label='Self-training')
plt.plot(x_values, label_propagation_score, label='Label Propagation')
plt.legend()
plt.ylabel("Accuracy")
plt.title("Comparision of classifiers on limited labeled data")
plt.xlabel("Amount of Labeled Data")
plt.show()
