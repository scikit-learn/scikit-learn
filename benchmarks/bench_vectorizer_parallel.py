# Author: Aldrian Obaja <aldrian_math@yahoo.co.id>
# License: BSD 3 clause
import time

import os
import pylab as pl
import nltk

from sklearn.feature_extraction.text import CountVectorizer
from sklearn.utils import check_random_state

import sys
class Unbuffered():
    def __init__(self, stream):
        self.stream = stream

    def write(self, data):
        self.stream.write(data)
        self.stream.flush()

    def __getattr__(self, attr):
        return getattr(self.stream, attr)
sys.stdout = Unbuffered(sys.stdout)

stemmer = nltk.stem.porter.PorterStemmer()

def filename_analyzer(filename):
    with open(filename,'r') as infile:
        return [stemmer.stem(word.lower()) for word in nltk.wordpunct_tokenize(infile.read())]

def string_analyzer(content):
    return [stemmer.stem(word.lower()) for word in nltk.wordpunct_tokenize(content)]

def plot(Vectorizer, corpora, analyzers):
    one_core = []
    multi_core = []

    for corpus, analyzer in zip(corpora, analyzers):
        vectorizer = Vectorizer(analyzer=analyzer)
        print 'Testing one-core...',
        start = time.time()
        vectorizer.fit_transform(corpus)
        one_core.append(time.time() - start)
        print 'Done in %.3fs' % one_core[-1]

        vectorizer = Vectorizer(analyzer=analyzer, n_jobs=-1)
        print 'Testing multi-core...',
        start = time.time()
        vectorizer.fit_transform(corpus)
        multi_core.append(time.time() - start)
        print 'Done in %.3fs' % multi_core[-1]

    pl.figure('scikit-learn parallel %s benchmark results' % Vectorizer.__name__)
    pl.plot([0,1], one_core, label="one core")
    pl.plot([0,1], multi_core, label="multi core")
    pl.xlabel('With/Without file loading')
    pl.ylabel('Time (s)')
    pl.title('Parallel %s' % Vectorizer.__name__)
    pl.legend()

list_of_filenames = [os.path.join(nltk.corpus.brown.root,fileid) for fileid in nltk.corpus.brown.fileids()]
list_of_contents = []
for filename in list_of_filenames:
    with open(filename,'r') as infile:
        list_of_contents.append(infile.read())
plot(CountVectorizer, [list_of_filenames, list_of_contents], [filename_analyzer, string_analyzer])

pl.show()
