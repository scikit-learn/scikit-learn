# Author: Aldrian Obaja <aldrian_math@yahoo.co.id>
# License: BSD 3 clause
#
# Benchmark the performance of Parallel version of CountVectorizer
#
# This will run the test on Brown corpus, containing 500 docs and
# 1,161,192 words.
#
# The result should be something like the following:
#
# $ python bench_vectorizer_parallel.py
# Testing one-core... Done in 38.760s
# Testing multi-core... Done in 22.352s
# Testing one-core... Done in 22.559s
# Testing multi-core... Done in 37.539s
# Testing one-core... Done in 0.009s
# Testing multi-core... Done in 0.138s
# Testing one-core... Done in 0.025s
# Testing multi-core... Done in 0.141s
#
# The first two do the benchmark on Brown corpus with heavy analyzer,
# which include stemming.
#
# The next two do the benchmark on Brown corpus repeated 7 times with
# light analyzer, which is the default when using analyzer='word'
#
# The next two do the benchmark on small dataset containing 35 docs
#
# The last two do the benchmark on small dataset containing 101 docs
import time

import os
import pylab as pl
import nltk

from sklearn.feature_extraction.text import CountVectorizer

import sys


class Unbuffered(object):

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
    with open(filename, 'r') as infile:
        return [stemmer.stem(word.lower())
                for word in nltk.wordpunct_tokenize(infile.read())]


def string_analyzer(content):
    return [stemmer.stem(word.lower())
            for word in nltk.wordpunct_tokenize(content)]


def plot(Vectorizer, corpora, analyzers,
         title='Parallel test at %s' % time.time(),
         xlabel='X-axis'):
    one_core = []
    multi_core = []

    for corpus, analyzer in zip(corpora, analyzers):
        vectorizer = Vectorizer(analyzer=analyzer)
        print 'Testing one-core...',
        start = time.time()
        vectorizer.fit_transform(corpus)
        one_core.append(time.time() - start)
        print 'Done in %.3fs' % one_core[-1]

        vectorizer = Vectorizer(analyzer=analyzer, n_jobs=2)
        print 'Testing multi-core...',
        start = time.time()
        vectorizer.fit_transform(corpus)
        multi_core.append(time.time() - start)
        print 'Done in %.3fs' % multi_core[-1]

    pl.figure(title)
    pl.plot(range(len(corpora)), one_core, label="one core")
    pl.plot(range(len(corpora)), multi_core, label="multi core")
    pl.xlabel(xlabel)
    pl.ylabel('Time (s)')
    pl.title('Parallel %s' % Vectorizer.__name__)
    pl.legend()

list_of_filenames = [os.path.join(nltk.corpus.brown.root, fileid)
                     for fileid in nltk.corpus.brown.fileids()]
list_of_contents = []
for filename in list_of_filenames:
    with open(filename, 'r') as infile:
        list_of_contents.append(infile.read())
plot(CountVectorizer,
     [list_of_filenames, list_of_contents * 7],
     [filename_analyzer, 'word'],
     title='CountVectorizer parallel benchmark on heavy/light analyzers',
     xlabel='Heavy Analyzer (500 docs) - Light Analyzer (3500 docs)')

plot(CountVectorizer,
     [['This is a very small document'] * 35,
         ['This is another small set of document'] * 101],
     [string_analyzer, string_analyzer],
     title='CountVectorizer parallel benchmark on small dataset',
     xlabel='35 small documents - 101 small documents')

pl.show()
