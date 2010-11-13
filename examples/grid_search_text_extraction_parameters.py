"""
==========================================================
Sample pipeline for text feature extraction and evaluation
==========================================================

The dataset used in this example is the 20 newsgroups dataset which will be
automatically downloaded and then cached and reused for the document
classification example.

You can adjust the number of categories by giving there name to the dataset
loader or setting them to None to get the 20 of them.
"""
print __doc__

# Author: Olivier Grisel <olivier.grisel@ensta.org>
#         Peter Prettenhofer <peter.prettenhofer@gmail.com>
#         Mathieu Blondel <mathieu@mblondel.org>
# License: Simplified BSD

from time import time
import os

from scikits.learn.datasets import load_files
from scikits.learn.feature_extraction.text.sparse import CountVectorizer
from scikits.learn.feature_extraction.text.sparse import TfidfTransformer
from scikits.learn.sgd.sparse import SGD
from scikits.learn.grid_search import GridSearchCV
from scikits.learn.pipeline import Pipeline

################################################################################
# Download the data, if not already on disk
url = "http://people.csail.mit.edu/jrennie/20Newsgroups/20news-18828.tar.gz"
archive_name = "20news-18828.tar.gz"

if not os.path.exists(archive_name[:-7]):
    if not os.path.exists(archive_name):
        import urllib
        print "Downloading data, please Wait (14MB)..."
        print url
        opener = urllib.urlopen(url)
        open(archive_name, 'wb').write(opener.read())
        print

    import tarfile
    print "Decompressiong the archive: " + archive_name
    tarfile.open(archive_name, "r:gz").extractall()
    print


################################################################################
# Load some categories from the training set
categories = [
    'alt.atheism',
    'talk.religion.misc',
]
# Uncomment the following to do the analysis on all the categories
#categories = None

print "Loading 20 newsgroups dataset for categories:"
print categories

data = load_files('20news-18828', categories=categories)
print "%d documents" % len(data.filenames)
print "%d categories" % len(data.target_names)
print

################################################################################
# define a pipeline combining a text feature extractor with a simple
# classifier
pipeline = Pipeline([
    ('vect', CountVectorizer()),
    ('tfidf', TfidfTransformer()),
    ('clf', SGD()),
])

parameters = {
#    'vect__max_df': (1.0, 0.7),
#    'vect__max_features': (10000, None),
    'vect__analyzer__max_n': (1, 2),
    'clf__alpha': (0.0001, 0.00001),
    'clf__penalty': ('l2', 'elasticnet'),
}

# find the best parameters for both the feature extraction and the
# classifier
grid_search = GridSearchCV(pipeline, parameters, n_jobs=-1)

# cross-validation doesn't work if the length of the data is not known,
# hence use lists instead of iterators
text_docs = [file(f).read() for f in data.filenames]

print "Performing grid search..."
print "pipeline:", [name for name, _ in pipeline.steps]
print "parameters:"
print parameters
t0 = time()
grid_search.fit(text_docs, data.target)
print "done in %0.3fs" % (time() - t0)
print

print "Best score: %0.3f" % grid_search.best_score
print "Best parameters set:"
best_parameters = grid_search.best_estimator._get_params()
for param_name in sorted(parameters.keys()):
    print "%s: %r" % (param_name, best_parameters[param_name])

