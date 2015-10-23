from sklearn import datasets
from sklearn.cross_validation import train_test_split
from sklearn.grid_search import AutomaticGridSearchCV
from sklearn.svm import SVC
from sklearn.linear_model.logistic import LogisticRegression

print(__doc__)

# Loading the Digits dataset
digits = datasets.load_digits()

# To apply a classifier on this data, we need to flatten the image, to
# turn the data in a (samples, feature) matrix:
n_samples = len(digits.images)
X = digits.images.reshape((n_samples, -1))
y = digits.target

# Split the dataset in two equal parts
X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.5, random_state=0)

print("# Tuning hyper-parameters\n")

classifiers = [SVC(), LogisticRegression()]

for clf in classifiers:
    clf_grid = AutomaticGridSearchCV(clf, cv=5, scoring='accuracy')
    clf_grid.fit(X_train, y_train)
    
    print "Best score for %s" % clf.__class__.__name__
    print clf_grid.best_score_
