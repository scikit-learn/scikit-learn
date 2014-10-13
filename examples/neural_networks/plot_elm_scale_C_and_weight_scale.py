"""
=================================================================
Impact of ELMs hyperparameters C and weight_scale on toy datasets
=================================================================

This illustrates how varying the regularization terms weight_scale, and C can
affect the regularity or smoothness of an ELM's decision function.

This example generates two plots, corresponding to each of the hyperparameters.
They control regularization in that they constrain the coefficients of the
decision function.

If, for example, there is high bias - as can be indicated by a high training
error, increasing C or weight_scale would decrease bias and therefore reduce
underfitting. But it is important not to increase their values so much that
they lead to high variance which indicates overfitting. Therefore, it is best
to find a balance between bias and variance (or a value that is not too high
nor too low) by testing a range of values as seen in this example.

"""
print(__doc__)


# Author: Issam H. Laradji
# License: BSD 3 clause

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from sklearn.cross_validation import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import make_moons, make_circles, make_classification
from sklearn.neural_network import ELMClassifier

h = .02  # step size in the mesh
rng = np.random.RandomState(1)

parameters = {}

parameters['C'] = np.logspace(-4, 4, 5)
parameters['weight_scale'] = np.logspace(-2, 2, 5)

default_parameters = ELMClassifier().get_params()

for param_name, param_values in sorted(parameters.items()):

    names = []
    classifiers = []
    for param_value in param_values:
        names.append("%s=%s" % (param_name, param_value))
        clf = ELMClassifier(n_hidden=100, random_state=rng)
        setattr(clf, param_name, param_value)
        classifiers.append(clf)

    X, y = make_classification(n_features=2, n_redundant=0, n_informative=2,
                               random_state=rng,
                               n_clusters_per_class=1)

    X += 2 * rng.uniform(size=X.shape)
    linearly_separable = (X, y)

    datasets = [make_moons(noise=0.3, random_state=rng),
                make_circles(noise=0.2, factor=0.5, random_state=rng),
                linearly_separable]

    figure = plt.figure(figsize=(17, 9))

    i = 1
    # iterate over datasets
    for X, y in datasets:
        # initialize standard scaler
        scaler = StandardScaler()

        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.4,
                                                            random_state=1)
        # Compute the mean and standard deviation of each feature of the
        # training set and scale the training set
        X_train = scaler.fit_transform(X_train)

        # Using the same mean and standard deviation, scale the testing set
        X_test = scaler.transform(X_test)

        x_min, x_max = X[:, 0].min() - .5, X[:, 0].max() + .5
        y_min, y_max = X[:, 1].min() - .5, X[:, 1].max() + .5
        xx, yy = np.meshgrid(np.arange(x_min, x_max, h),
                             np.arange(y_min, y_max, h))

        # just plot the dataset first
        cm_bright = ListedColormap(['#FF0000', '#0000FF'])
        ax = plt.subplot(len(datasets), len(classifiers) + 1, i)
        # Plot the training points
        ax.scatter(X_train[:, 0], X_train[:, 1], c=y_train, cmap=cm_bright)
        # and testing points
        ax.scatter(X_test[:, 0], X_test[:, 1], c=y_test, cmap=cm_bright,
                   alpha=0.6)
        ax.set_xlim(xx.min(), xx.max())
        ax.set_ylim(yy.min(), yy.max())
        ax.set_xticks(())
        ax.set_yticks(())
        i += 1

        # iterate over classifiers
        for name, clf in zip(names, classifiers):
            ax = plt.subplot(len(datasets), len(classifiers) + 1, i)
            clf.fit(X_train, y_train)
            score = clf.score(X_test, y_test)

            # Plot the decision boundary.
            Z = clf.predict_proba(np.c_[xx.ravel(), yy.ravel()])[:, 1]

            # Put the result into a color plot
            Z = Z.reshape(xx.shape)
            ax.contourf(xx, yy, Z, cmap=plt.cm.RdBu, alpha=.8)

            # Plot also the training points
            ax.scatter(X_train[:, 0], X_train[:, 1], c=y_train, cmap=cm_bright)
            # and testing points
            ax.scatter(X_test[:, 0], X_test[:, 1], c=y_test, cmap=cm_bright,
                       alpha=0.6)

            ax.set_xlim(xx.min(), xx.max())
            ax.set_ylim(yy.min(), yy.max())
            ax.set_xticks(())
            ax.set_yticks(())
            ax.set_title(name)
            ax.text(xx.max() - .3, yy.min() + .3, ('%.2f' % score).lstrip('0'),
                    size=15, horizontalalignment='right')
            i += 1
    figure.subplots_adjust(left=.02, right=.98)

    other_param_name = 'C' if param_name == 'weight_scale' else 'weight_scale'
    other_param_value = default_parameters[other_param_name]
    figure.suptitle("Impact of varying {param_name} for fixed"
                    " {other_param_name}={other_param_value}".format(
                        param_name=param_name,
                        other_param_name=other_param_name,
                        other_param_value=other_param_value))
plt.show()
