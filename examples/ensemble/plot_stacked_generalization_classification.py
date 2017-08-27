"""
==========================================================
Stacked generalization applied to a classification problem
==========================================================

This example shows how stacked generalization can be used to combine several
classifiers into a single stacked one. We used only two of the features to be
able to visualize it better. Because of that, the model performs badly, but
what's interesting here is to see how the decision boundaries are combined by
the stacking framework.
"""

# utils
import numpy as np
import matplotlib.pyplot as plt

# import base regressors
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier

# stacking api
from sklearn.ensemble import make_stack_layer
from sklearn.pipeline import Pipeline

# dataset
from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split

RANDOM_SEED = 89555

# setup plot
SUBPLOT_X = 1
SUBPLOT_Y = 4
SUBPLOT_OFFSET = (SUBPLOT_X * 10 + SUBPLOT_Y) * 10
EXAMPLE_INDEX = 1

plt.figure(EXAMPLE_INDEX, figsize=(9, 3))

# prepare data
X, y = load_iris(return_X_y=True)
X = X[:, :2]  # use only two features to make it easier to visualize
X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    random_state=RANDOM_SEED)

# Base classifier
svc = SVC(C=1e3, gamma=1e-2, probability=True, random_state=RANDOM_SEED)
et = ExtraTreesClassifier(max_depth=None, n_estimators=30,
                          min_samples_split=2, random_state=RANDOM_SEED)
rf = RandomForestClassifier(max_depth=None, n_estimators=30,
                            random_state=RANDOM_SEED)

base_classifiers = [("SVC", svc),
                    ("ExtraTrees", et),
                    ("RandomForest", rf)]


def plot_model(name, model, plot_idx):
    model.fit(X_train, y_train)

    plt.subplot(SUBPLOT_OFFSET + plot_idx)
    plt.title(name)
    x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
    y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
    xx, yy = np.meshgrid(np.arange(x_min, x_max, .2),
                         np.arange(y_min, y_max, .2))
    z = model.predict(np.c_[xx.ravel(), yy.ravel()]).reshape(xx.shape)
    plt.contourf(xx, yy, z, cmap=plt.cm.RdYlBu, alpha=.5)
    plt.scatter(X[:, 0], X[:, 1], c=y, s=15, cmap=plt.cm.RdYlBu)
    plt.xticks(())
    plt.yticks(())


for i, (name, classifier) in enumerate(base_classifiers):
    plot_model(name, classifier, 1 + i)

# Stacked ensemble: we use the base classifiers as features for another model


layer0 = make_stack_layer(base_classifiers)

X2_train = layer0.fit_transform(X_train, y_train)

combiner = LogisticRegression(random_state=RANDOM_SEED)
final_classifier = Pipeline([('layer0', layer0),
                            ('layer1', combiner)])

plot_model('Stacked classifiers', final_classifier, SUBPLOT_Y)

plt.tight_layout()
plt.show()
