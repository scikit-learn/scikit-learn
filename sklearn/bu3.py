import numpy as np
import matplotlib.pyplot as plt

from sklearn import svm, datasets
from sklearn.model_selection import train_test_split
from sklearn.metrics import plot_confusion_matrix, confusion_matrix, ConfusionMatrixDisplay
from sklearn.utils.multiclass import unique_labels

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

    y_pred = classifier.predict(X_test)

    display_labels = unique_labels(y_test, y_pred)

    cm = confusion_matrix(y_test, y_pred,
                          labels=display_labels, normalize=normalize)

    disp = ConfusionMatrixDisplay(confusion_matrix=cm,
                                  display_labels=display_labels, font_size=20)

    # optional way
    #disp.font_size = 40

    disp.plot(cmap=plt.cm.Blues)

    disp.ax_.set_title(title)

plt.show()