from scikits.learn.cross_val import StratifiedKFold
from scikits.learn.datasets import fetch_lfw_people
from scikits.learn.grid_search import GridSearchCV
from scikits.learn.decomposition import RandomizedPCA
from scikits.learn.svm import SVC

# Download the data, if not already on disk and load it as numpy arrays
lfw_people = fetch_lfw_people(min_faces_per_person=70, resize=0.4)
# reshape the data using the traditional (n_samples, n_features) shape
faces = lfw_people.data
n_samples, h, w = faces.shape
X = faces.reshape((n_samples, h * w))
n_features = X.shape[1]

# the label to predict is the id of the person
y = lfw_people.target
target_names = lfw_people.target_names

# split into a training and testing set
train, test = iter(StratifiedKFold(y, k=4)).next()
X_train, X_test = X[train], X[test]
y_train, y_test = y[train], y[test]

# Compute a PCA (eigenfaces) on the face dataset (treated as unlabeled
# dataset): unsupervised feature extraction / dimensionality reduction
n_components = 150
pca = RandomizedPCA(n_components=n_components, whiten=True).fit(X_train)
eigenfaces = pca.components_.reshape((n_components, h, w))

X_train_pca = pca.transform(X_train)
X_test_pca = pca.transform(X_test)

# Train a SVM classification model
param_grid = dict(C=[1, 5, 10, 50, 100],
                  gamma=[0.0001, 0.0005, 0.001, 0.005, 0.01, 0.1])
clf = GridSearchCV(SVC(kernel='rbf'), param_grid,
                   fit_params={'class_weight': 'auto'},
                   verbose=1)
clf = clf.fit(X_train_pca, y_train)
print clf.best_estimator

# Quantitative evaluation of the model quality on the test set
from scikits.learn import metrics
y_pred = clf.predict(X_test_pca)
print metrics.classification_report(y_test, y_pred, target_names=target_names)
print metrics.confusion_matrix(y_test, y_pred,
                               labels=range(len(target_names)))


# Plot the results
import pylab as pl
for index, (img, label_true, label_pred) in enumerate(
                zip(X_test[:8], y_test[:8], y_pred[:8])):
    pl.subplot(2, 4, index+1).imshow(img.reshape(h, w), cmap=pl.cm.gray)
    pl.title('%s, prediction: %s' % (label_true, label_pred))

