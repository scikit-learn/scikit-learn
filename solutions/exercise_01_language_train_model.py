"""Build a language detector model"""
# Author: Olivier Grisel <olivier.grisel@ensta.org>
# License: Simplified BSD

import sys

from scikits.learn.feature_extraction.text.sparse import Vectorizer
from scikits.learn.feature_extraction.text import CharNGramAnalyzer
from scikits.learn.svm.sparse import LinearSVC
from scikits.learn.datasets import load_files
from scikits.learn import metrics

#
# New preprocessor better suited for language id than the default
# preprocessor
#

class LowerCasePreprocessor(object):
    """Simple preprocessor that should be available by default"""

    def preprocess(self, unicode_content):
        return unicode_content.lower()

    def __repr__(self):
        return "LowerCasePreprocessor()"

#
# The real code starts here
#


# the training data folder must be passed as first argument
languages_data_folder = sys.argv[1]
dataset = load_files(languages_data_folder)

# split the dataset in training and test set:
n_samples_total = dataset.filenames.shape[0]

filenames_train = dataset.filenames[:n_samples_total/2]
filenames_test = dataset.filenames[n_samples_total/2:]

y_train = dataset.target[:n_samples_total/2]
y_test = dataset.target[n_samples_total/2:]


# Build a an analyzer that split strings into sequence of 1 to 3 characters
# after using the previous preprocessor
analyzer = CharNGramAnalyzer(
    min_n=1,
    max_n=3,
    preprocessor=LowerCasePreprocessor(),
)

# Build a vectorizer using the analyzer, learn the mapping from feature name to
# feature id on the training data while transforming it. The use the fitted
# vectorizer on the test data
vectorizer = Vectorizer(analyzer=analyzer, use_idf=False)

X_train = vectorizer.fit_transform((open(f) for f in filenames_train))
X_test = vectorizer.transform((open(f) for f in filenames_test))

# Build a linear classifier and train it on the training set
clf = LinearSVC(loss='l2', penalty='l1', dual=False, C=100)
clf.fit(X_train, y_train)

# Predict the outcome on the testing set
y_predicted = clf.predict(X_test)

# Print the classification report
print metrics.classification_report(y_test, y_predicted,
                                    class_names=dataset.target_names)

# Plot the confusion matrix
cm = metrics.confusion_matrix(y_test, y_predicted)
print cm

# import pylab as pl
#pl.matshow(cm)
#pl.show()

# Predict the result on some short new sentences:
sentences = [
    u'This is a language detection test.',
    u'Ceci est un test de d\xe9tection de la langue.',
    u'Dies ist ein Test, um die Sprache zu erkennen.',
]
vectors = vectorizer.transform(sentences)
predicted = clf.predict(vectors)

for s, p in zip(sentences, predicted):
    print u'The language of "%s" is "%s"' % (s, dataset.target_names[p])

