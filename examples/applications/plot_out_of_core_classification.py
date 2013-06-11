"""
======================================================
Out-of-core classification of text documents
======================================================

This is an example showing how scikit-learn can be used for classification
using an out-of-core approach. This example uses a `HashingVectorizer`
and a classifier supporting `partial_fit ` to limit memory consumption.

The dataset used in this example is Reuters-21578 as provided by the UCI ML
repository. It will be automatically downloaded and uncompressed in the current
directory on first run.

The plot represents the evolution of classification accuracy with the number
of mini-batches fed to the classifier.
"""

# Author: Eustache Diemert <eustache@diemert.fr>
# License: BSD 3 clause

from __future__ import print_function

import time
import random
import re
import os.path
import fnmatch
import sgmllib
import urllib
import tarfile

import numpy as np
import pylab as pl

from sklearn.feature_extraction.text import HashingVectorizer
from sklearn.linear_model.stochastic_gradient import SGDClassifier

###############################################################################
# Reuters Dataset related routines
###############################################################################


class ReutersParser(sgmllib.SGMLParser):
    """Utility class to parse a SGML file and yield documents one at a time."""
    def __init__(self, verbose=0):
        sgmllib.SGMLParser.__init__(self, verbose)
        self._reset()

    def _reset(self):
        self.in_title = 0
        self.in_body = 0
        self.in_topics = 0
        self.in_topic_d = 0
        self.title = ""
        self.body = ""
        self.topics = []
        self.topic_d = ""

    def parse(self, fd):
        self.docs = []
        for chunk in fd:
            self.feed(chunk)
            for doc in self.docs:
                yield doc
            self.docs = []
        self.close()

    def handle_data(self, data):
        if self.in_body:
            self.body += data
        elif self.in_title:
            self.title += data
        elif self.in_topic_d:
            self.topic_d += data

    def start_reuters(self, attributes):
        pass

    def end_reuters(self):
        self.body = re.sub(r'\s+', r' ', self.body)
        self.docs.append({'title': self.title,
                          'body': self.body,
                          'topics': self.topics})
        self._reset()

    def start_title(self, attributes):
        self.in_title = 1

    def end_title(self):
        self.in_title = 0

    def start_body(self, attributes):
        self.in_body = 1

    def end_body(self):
        self.in_body = 0

    def start_topics(self, attributes):
        self.in_topics = 1

    def end_topics(self):
        self.in_topics = 0

    def start_d(self, attributes):
        self.in_topic_d = 1

    def end_d(self):
        self.in_topic_d = 0
        self.topics.append(self.topic_d)
        self.topic_d = ""


class ReutersStreamReader():
    """Iterate over documents stored in a local directory holding the Reuters
    dataset.

    The Reuters archive will automatically be downloaded and uncompressed if
    the `data_path` directory does not exist.

    Documents are represented as dictionaries with 'body' (str),
    'title' (str), 'topics' (list(str)) keys.
    """

    DOWNLOAD_URL = ('http://archive.ics.uci.edu/ml/machine-learning-databases/'
                    'reuters21578-mld/reuters21578.tar.gz')
    ARCHIVE_FILENAME = 'reuters21578.tar.gz'

    def __init__(self, data_path):
        self.data_path = data_path
        if not os.path.exists(self.data_path):
            self.download_dataset()

    def download_dataset(self):
        print("downloading dataset (once and for all) into %s" %
              self.data_path)
        os.mkdir(self.data_path)

        def progress(blocknum, bs, size):
            total_sz_mb = '%.2f MB' % (size / 1e6)
            current_sz_mb = '%.2f MB' % ((blocknum * bs) / 1e6)
            print('\rdownloaded %s / %s' % (current_sz_mb, total_sz_mb),
                  end='')
        urllib.urlretrieve(self.DOWNLOAD_URL,
                           filename=os.path.join(self.data_path,
                                                 self.ARCHIVE_FILENAME),
                           reporthook=progress)
        print('\r', end='')
        print("untaring data ...")
        tfile = tarfile.open(os.path.join(self.data_path,
                                          self.ARCHIVE_FILENAME),
                             'r:gz')
        tfile.extractall(self.data_path)
        print("done !")

    def iterdocs(self):
        for root, _dirnames, filenames in os.walk(self.data_path):
            for filename in fnmatch.filter(filenames, '*.sgm'):
                path = os.path.join(root, filename)
                parser = ReutersParser()
                for doc in parser.parse(open(path)):
                    yield doc


###############################################################################
# Main
###############################################################################
"""Create the hasher and limit the nber of features to a reasonable maximum."""
hasher = HashingVectorizer(charset_error='ignore', n_features=2 ** 18)

"""Create an online classifier i.e. supporting `partial_fit()`."""
classifier = SGDClassifier()

"""Create the data_streamer that parses Reuters SGML files and iterates on
documents as a stream."""
data_streamer = ReutersStreamReader('./reuters/')

"""Here we propose to learn a binary classification between the positive class
and all other documents."""
all_classes = np.array([0, 1])
positive_class = 'acq'

"""We will feed the classifier with mini-batches of 100 documents; this means
we have at most 100 docs in memory at any time."""
chunk = []
chunk_sz = 100

stats = {'n_train': 0, 'n_test': 0, 'n_train_pos': 0, 'n_test_pos': 0,
         'accuracy': 0.0, 'accuracy_history': [(0, 0)], 't0': time.time(),
         'runtime_history': [(0, 0)]}


def progress(stats):
    """Reports progress information."""
    s = "%(n_train)d train docs (%(n_train_pos)d positive) " % stats
    s += "%(n_test)d test docs (%(n_test_pos)d positive) " % stats
    s += "accuracy: %(accuracy)f " % stats
    s += "in %.2fs" % (time.time() - stats['t0'])
    return s

# Main loop : iterate over documents read by the streamer
for i, doc in enumerate(data_streamer.iterdocs()):

    if i and not i % 10:
        # Print progress information
        print("\r%s" % progress(stats), end='')

    # Discard invalid documents
    if not len(doc['topics']):
        continue

    # Read documents until chunk full
    if len(chunk) < chunk_sz:
        classid = int(positive_class in doc['topics'])
        chunk.append((doc['title'] + '\n\n' + doc['body'], classid))
        continue

    # When chunk is full, create data matrix using the HashingVectorizer
    documents, topics = zip(*chunk)
    y = np.array(topics)
    X = hasher.transform(documents)
    chunk = []

    # Once every 10 chunks or so, test accuracy.
    if random.random() < 0.1:
        stats['n_test'] += len(documents)
        stats['n_test_pos'] += sum(topics)
        stats['accuracy'] = classifier.score(X, y)
        stats['accuracy_history'].append((stats['accuracy'], stats['n_train']))
        stats['runtime_history'].append((stats['accuracy'],
                                         time.time() - stats['t0']))
        continue

    # Learn from the current chunk.
    stats['n_train'] += len(documents)
    stats['n_train_pos'] += sum(topics)
    classifier.partial_fit(X, y, classes=all_classes)

print()

# Plot accuracy evolution with #examples
pl.figure(1)
pl.subplots_adjust(hspace=0.5)
pl.subplot(211)
pl.title('Classification accuracy as a function of #examples seen')
pl.xlabel('# training examples')
pl.ylabel('Accuracy')
y, x = zip(*stats['accuracy_history'])
x = np.array(x)
y = np.array(y)
pl.grid(True)
pl.plot(x, y)
# Plot accuracy evolution with runtime
pl.subplot(212)
pl.title('Classification accuracy as a function of runtime')
pl.xlabel('Runtime (s)')
pl.ylabel('Accuracy')
y, x = zip(*stats['runtime_history'])
x = np.array(x)
y = np.array(y)
pl.grid(True)
pl.plot(x, y)
pl.show()
