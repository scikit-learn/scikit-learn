"""
======================================================
Out-of-core classification of text documents 
======================================================

This is an example showing how scikit-learn can be used for classification 
using an out-of-core approach. This example uses a `HashingVectorizer`
and a classifier supporting `partial_fit ` to limit memory consumption. 

The dataset used in this example is Reuters-21578 as provided by the UCI ML 
repository. It should be downloaded and uncompressed in the current directory.
e.g. wget http://archive.ics.uci.edu/ml/machine-learning-databases/reuters21578-mld/reuters21578.tar.gz && tar xvzf reuters21578.tar.gz

The plot represents the evolution of classification accuracy with the number
of mini-batches fed to the classifier.
"""

# Author: Eustache Diemert <eustache@diemert.fr>
# License: BSD 3 clause

import sys
import time
import random
import re
import os.path
import fnmatch
import sgmllib
from collections import defaultdict

import numpy as np
import pylab as pl

from sklearn.feature_extraction import FeatureHasher
from sklearn.linear_model.stochastic_gradient import SGDClassifier

################################################################################
# Reuters Dataset related routines
################################################################################

class ReutersTopics():
    """Utility class to read official topic names from the relevant file."""
    TOPICS_FILENAME = 'all-topics-strings.lc.txt'
    def __init__(self, topics_path):
        self.topics_ = open(topics_path).read().split('\n')
        self.topics_ = dict([ (self.topics_[i],i) for i in range(len(self.topics_)) ])

    def topic_ids(self):
        return self.topics_.values()



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
        self.docs.append({'title':self.title,
                          'body':self.body,
                          'topics':self.topics})
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
    """Reads documents from the directory where the Reuters dataset has been 
    uncompressed; will download the dataset if needed.
    
    Documents are represented as dictionaries with 'body' (str), 
    'title' (str), 'topics' (list(str)) keys. 
    """
    
    DOWNLOAD_URL = 'http://archive.ics.uci.edu/ml/machine-learning-databases/reuters21578-mld/reuters21578.tar.gz'
    ARCHIVE_FILENAME = 'reuters21578.tar.gz'
    
    def __init__(self, data_path):
        self.data_path = data_path
        if not os.path.exists(self.data_path):
            self.download_dataset()
        self.topics = ReutersTopics(os.path.join(data_path,
                                                 ReutersTopics.TOPICS_FILENAME))
        self.classes = self.topics.topic_ids()

    def download_dataset(self):
        print "downloading dataset (once and for all) into %s"%self.data_path
        os.mkdir(self.data_path)
        import urllib
        def progress(blocknum, bs, size):
            total_sz_mb = '%.2f MB'%(size/1e6)
            current_sz_mb = '%.2f MB'%((blocknum*bs)/1e6)
            print >>sys.stderr,'\rdownloaded %s / %s'%(current_sz_mb,
                                                       total_sz_mb),
        urllib.urlretrieve(self.DOWNLOAD_URL, 
                           filename=os.path.join(self.data_path,
                                                 self.ARCHIVE_FILENAME),
                           reporthook=progress)
        print >>sys.stderr, '\r'
        import tarfile
        print "untaring data ...",
        tfile = tarfile.open(os.path.join(self.data_path, 
                                          self.ARCHIVE_FILENAME), 
                             'r:gz')
        tfile.extractall(self.data_path)
        print "done !"

    def iterdocs(self):
        for root, _dirnames, filenames in os.walk(self.data_path):
            for filename in fnmatch.filter(filenames, '*.sgm'):
                path = os.path.join(root, filename)
                parser = ReutersParser()
                for doc in parser.parse(open(path)):
                    yield doc



################################################################################
# Feature extraction routines
################################################################################


def tokens(doc):
    """Extract tokens from doc.

    This uses a simple regex to break strings into tokens. For a more
    principled approach, see CountVectorizer or TfidfVectorizer.
    """
    return (tok.lower() for tok in re.findall(r"\w+", doc))

def token_freqs(doc, freq=None):
    """Extract a dict mapping tokens from doc to their frequencies."""
    if not freq:
        freq = defaultdict(int)
    for tok in tokens(doc):
        freq[tok] += 1
    return freq



################################################################################
# Main
################################################################################


"""Create the hasher and limit the nber of features to a reasonable maximum."""
hasher = FeatureHasher(n_features=2**18)

"""Create an online classifier i.e. supporting `partial_fit()`."""
classifier = SGDClassifier()

"""Create the data_streamer that parses Reuters SGML files and iterates on 
documents as a stream."""
data_streamer = ReutersStreamReader('./reuters/')

"""Here we propose to learn a binary classification between the positive class
and all other documents."""
all_classes = np.array([0,1])
positive_class = 'acq'

"""We will feed the classifier with mini-batches of 100 documents; this means
we have at most 100 docs in memory at any time."""
chunk = []
chunk_sz = 100

stats = {'n_train':0,'n_test':0,'n_train_pos':0,'n_test_pos':0,'accuracy':0.0,
         'accuracy_history':[(0,0)],'t0':time.time()}

def progress(stats):
    """Reports progress information."""
    s = "%(n_train)d train docs (%(n_train_pos)d positive) "%stats
    s+= "%(n_test)d test docs (%(n_test_pos)d positive) "%stats
    s+= "accuracy: %(accuracy)f "%stats
    s+= "in %.2fs"%(time.time()-stats['t0'])
    return s

"""Main loop : iterate over documents read by the streamer."""
for i, doc in enumerate(data_streamer.iterdocs()):

    if i and not i % 10:
        """Print progress information."""
        print >>sys.stderr, "\r%s"%progress(stats),

    """Discard invalid documents."""
    if not len(doc['topics']):
        continue

    """Read documents until chunk full."""
    if len(chunk) < chunk_sz:
        freq = token_freqs(doc['title'])
        freq = token_freqs(doc['body'], freq)
        classid = int(positive_class in doc['topics'])
        chunk.append((freq, classid))
        continue

    """When chunk is full, create data matrix using the HashingVectorizer"""
    freqs, topics = zip(*chunk)
    y = np.array(topics)
    X = hasher.transform(freqs)
    chunk = []

    """Once every 10 chunks or so, test accuracy."""
    if random.random() < 0.1:
        stats['n_test'] += len(freqs)
        stats['n_test_pos'] += sum(topics)
        stats['accuracy'] = classifier.score(X, y)
        stats['accuracy_history'].append((stats['accuracy'],
                                          stats['n_train']))
        continue

    """Learn from the current chunk."""
    stats['n_train'] += len(freqs)
    stats['n_train_pos'] += sum(topics)
    classifier.partial_fit(X, 
                           y,
                           classes=all_classes)

print >>sys.stderr

"""Plot accuracy evolution with time."""
pl.figure()
pl.title('Classification accuracy as a function of #examples seen')
pl.xlabel('# training examples')
pl.ylabel('Accuracy')
y,x = zip(*stats['accuracy_history'])
x = np.array(x)
y = np.array(y)
pl.plot(x,y)
pl.show()

