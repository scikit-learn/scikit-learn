"""
===================================================================
Recognizing hand-written digits using Fastfood kernel approximation
===================================================================

This shows how the Fastfood kernel approximation compares to a dual and primal
support vector classifier. It is based on the plot_digits_classification
example of scikit-learn. The idea behind Fastfood is to map the data into a
feature space (approximation) and then run a linear classifier on the mapped
data.


"""

print __doc__

# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org>
# Modified By: Felix Maximilan Möller
# License: Simplified BSD

# Standard scientific Python imports
import numpy as np
import pylab as pl

# Import datasets, classifiers and performance metrics
from sklearn import datasets, svm, metrics
from sklearn.kernel_approximation import Fastfood

# The digits dataset
digits = datasets.load_digits()

# The data that we are interested in is made of 8x8 images of digits,
# let's have a look at the first 3 images, stored in the `images`
# attribute of the dataset. If we were working from image files, we
# could load them using pylab.imread. For these images know which
# digit they represent: it is given in the 'target' of the dataset.
for index, (image, label) in enumerate(zip(digits.images, digits.target)[:4]):
    pl.subplot(2, 4, index + 1)
    pl.axis('off')
    pl.imshow(image, cmap=pl.cm.gray_r, interpolation='nearest')
    pl.title('Training: %i' % label)

# To apply an classifier on this data, we need to flatten the image, to
# turn the data in a (samples, feature) matrix:
n_samples = len(digits.images)
data = digits.images.reshape((n_samples, -1))
gamma = .001
sigma = np.sqrt(1 / (2 * gamma))
number_of_features_to_generate = 1000
train__idx = range(n_samples / 2)
test__idx = range(n_samples / 2, n_samples)

# map data into featurespace
rbf_transform = Fastfood(sigma=sigma, n_components=number_of_features_to_generate)
data_transformed_train = rbf_transform.fit(data[train__idx]).transform(data[train__idx])
data_transformed_test = rbf_transform.transform(data[test__idx])

# Create a classifier: a support vector classifier
classifier = svm.SVC(gamma=gamma)
linear_classifier = svm.LinearSVC()
linear_classifier_transformation = svm.LinearSVC()

# We learn the digits on the first half of the digits
classifier.fit(data[train__idx], digits.target[train__idx])
linear_classifier.fit(data[train__idx], digits.target[train__idx])

# Run the linear classifier on the mapped data.
linear_classifier_transformation.fit(data_transformed_train, digits.target[train__idx])

# Now predict the value of the digit on the second half:
expected = digits.target[test__idx]
predicted = classifier.predict(data[test__idx])
predicted_linear = linear_classifier.predict(data[test__idx])
predicted_linear_transformed = linear_classifier_transformation.predict(data_transformed_test)

print "Classification report for dual classifier %s:\n%s\n" % (
    classifier, metrics.classification_report(expected, predicted))
print "Classification report for primal linear classifier %s:\n%s\n" % (
    linear_classifier, metrics.classification_report(expected, predicted_linear))
print "Classification report for primal transformation classifier %s:\n%s\n" % (
    linear_classifier_transformation, metrics.classification_report(expected, predicted_linear_transformed))

print "Confusion matrix for dual classifier:\n%s" % metrics.confusion_matrix(expected, predicted)
print "Confusion matrix for primal linear classifier:\n%s" % metrics.confusion_matrix(expected, predicted_linear)
print "Confusion matrix for for primal transformation classifier:\n%s" % metrics.confusion_matrix(expected, predicted_linear_transformed)

# assert_almost_equal(metrics.classification_report(expected, predicted),
#                     metrics.classification_report(expected, predicted_linear_transformed),
#                     decimal=1)

for index, (image, prediction) in enumerate(zip(digits.images[test__idx], predicted)[:4]):
    pl.subplot(2, 4, index + 5)
    pl.axis('off')
    pl.imshow(image, cmap=pl.cm.gray_r, interpolation='nearest')
    pl.title('Prediction: %i' % prediction)

pl.show()
