"""Build a sentiment analysis / polarity model"""
# Author: Olivier Grisel <olivier.grisel@ensta.org>
# License: Simplified BSD

import sys
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.svm.sparse import LinearSVC
from sklearn.pipeline import Pipeline
from sklearn.grid_search import GridSearchCV
from sklearn.datasets import load_files
from sklearn import metrics

#
# The real code starts here
#


# the training data folder must be passed as first argument
movie_reviews_data_folder = sys.argv[1]
dataset = load_files(movie_reviews_data_folder, shuffle=True, random_state=42)

# split the dataset in training and test set:
n_samples_total = dataset.filenames.shape[0]

split = (n_samples_total * 3) / 4

docs_train = dataset.data[:split]
docs_test = dataset.data[split:]

y_train = dataset.target[:split]
y_test = dataset.target[split:]

# Build a vectorizer / classifier pipeline using the previous analyzer

# TODO

# Build a grid search to find out whether unigrams or bigrams are more
# useful

parameters = {
    'vect__analyzer__max_n': (1, 2),
}

# TODO

## Predict the outcome on the testing set
#y_predicted = clf.predict(docs_test)
#
## Print the classification report
#print metrics.classification_report(y_test, y_predicted,
#                                    target_names=dataset.target_names)
#
## Plot the confusion matrix
#cm = metrics.confusion_matrix(y_test, y_predicted)
#print cm
#
## import pylab as pl
##pl.matshow(cm)
##pl.show()

