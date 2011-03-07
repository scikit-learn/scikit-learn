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

# TODO: define variables 'filenames_train' and 'filenames_test'
# TODO: define variables 'y_train' and 'y_test'


# Build a an analyzer that split strings into sequence of 1 to 3 characters
# using the previous preprocessor

# TODO: define a variable named analyzer


# Build a vectorizer using the analyzer, learn the mapping from feature name to
# feature id on the training data and then transform it into feature vectors.
# Then use the fitted vectorizer on the test data

# TODO: define a variable named 'vectorizer'
# TODO: define a variable named 'X_train'
# TODO: define a variable named 'X_test'

# XXX: Don't forget to read the content of the text files before feeding it to
# the vectorizer

# Build a linear classifier and train it on the training set

# TODO: define a variable named 'clf'

# Predict the outcome on the testing set

# TODO: define a variable named 'y_predicted'


#
# Evaluation of the quality of the predictions: uncomment the following when all
# of the above as been implemented
#

## Print the classification report
#
#print metrics.classification_report(y_test, y_predicted,
#                                    class_names=dataset.target_names)
#
## Print the confusion matrix
#
#cm = metrics.confusion_matrix(y_test, y_predicted)
#print cm
#
# Predict the result on some short new sentences:
#sentences = [
#    u'This is a language detection test.',
#    u'Ceci est un test de d\xe9tection de la langue.',
#    u'Dies ist ein Test, um die Sprache zu erkennen.',
#]
#vectors = vectorizer.transform(sentences)
#predicted = clf.predict(vectors)
#
#for s, p in zip(sentences, predicted):
#    print u'The language of "%s" is "%s"' % (s, dataset.target_names[p])

