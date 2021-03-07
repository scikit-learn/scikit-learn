import numpy as np
import matplotlib.pyplot as plt

from sklearn import svm, datasets
from sklearn.model_selection import train_test_split
from sklearn.metrics import plot_confusion_matrix
#import confusion_matrix
from sklearn.metrics import ConfusionMatrixDisplay, confusion_matrix
from sklearn.utils.multiclass import unique_labels
from sklearn.datasets import make_classification
from sklearn.svm import SVC



# X, y = make_classification(random_state=0)
# X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
# clf = SVC(random_state=0)
# clf.fit(X_train, y_train)

# y_pred = clf.predict(X)
# display_labels = unique_labels(y_test, y_pred)
# cm = confusion_matrix(y_test, y_pred, sample_weight=None,labels=display_labels, normalize=True)
# ConfusionMatrixDisplay(cm, font_size=50).from_estimator(clf, X_test, y_test)

# plt.show()

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

print(X_test)

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
    print(display_labels)
    print(y_test)
    #y_test.put("versicolor")
    cm = confusion_matrix(y_test, y_pred, labels=None,sample_weight=None,
                          normalize=normalize)

    display_labels = unique_labels(y_test, y_pred)

    print(display_labels)
    disp = ConfusionMatrixDisplay(confusion_matrix=cm,
                                  display_labels=display_labels,val=40)
    
    #print(disp.font_size)
    print(disp.display_labels)
    disp.font_size = 40
    #print(disp.font_size)

    #ConfusionMatrixDisplay.from_estimator(classifier, X_test, y_test)
    disp.plot(cmap=plt.cm.Blues)
    #disp.ax.set_title(title)


plt.show()