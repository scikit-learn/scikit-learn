import numpy as np
import matplotlib.pyplot as plt

from sklearn import svm, datasets
from sklearn.model_selection import train_test_split
from sklearn.metrics import plot_confusion_matrix

# import some data to play with
iris = datasets.load_iris()
X = iris.data
y = iris.target
class_names = iris.target_names

# Split the data into a training set and a test set
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

# Run classifier, using a model that is too regularized (C too low) to see
# the impact on the results
classifier = svm.SVC(kernel='linear', C=0.01).fit(X_train, y_train)

np.set_printoptions(precision=2)

# Plot non-normalized confusion matrix
titles_options = [("Confusion matrix, without normalization", 'true'),
                  ("Normalized confusion matrix", 'true')]
for title, normalize in titles_options:

    display_labels = class_names
    include_values = True
    cmap = plt.cm.Blues
    if not cmap :
        cmap = 'viridis'

    y_pred = classifier.predict(X_test)
    cm = confusion_matrix(y_test, y_pred, sample_weight=None,
                          labels=display_labels, normalize=normalize)

    '''if display_labels is None:
        if labels is None:
            display_labels = unique_labels(y_test, y_pred)
        else:
            display_labels = labels'''

    disp = ConfusionMatrixDisplay(confusion_matrix=cm,
                                  display_labels=display_labels, font_size=30)

    disp.from_predictions(classifier, y_test, y_pred, display_labels=class_names,
                                 cmap=plt.cm.Blues,
                                 normalize=normalize)

    '''disp.plot(include_values=include_values,
                     cmap=plt.cm.Blues, ax=None, xticks_rotation='horizontal',
                     values_format=None, colorbar=True)'''

   ''' disp = plot_confusion_matrix(classifier, X_test, y_test,
                                 display_labels=class_names,
                                 cmap=plt.cm.Blues,
                                 normalize=normalize, font_size=None)'''
    disp.ax_.set_title(title)


plt.show()