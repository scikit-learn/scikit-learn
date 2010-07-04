# Author: Olivier Grisel <olivier.grisel@ensta.org>
#
# License: BSD Style.
"""Utilities to build feature vectors from text documents"""

import re
import unicodedata
import numpy as np


def strip_accents(s):
    """Transform accentuated unicode symbols into their simple counterpart"""
    return ''.join((c for c in unicodedata.normalize('NFD', s)
                    if unicodedata.category(c) != 'Mn'))


class SimpleAnalyzer(object):
    """Simple analyzer: transform a text document into a sequence of tokens

    This simple implementation does:
        - lower case conversion
        - unicode accents removal
        - token extraction using unicode regexp word bounderies for token of
          minimum size of 2 symbols
    """

    token_pattern = re.compile(r"\b\w\w+\b", re.U)

    def analyze(self, text_document):
        if isinstance(text_document, str):
            text_document = text_document.decode("utf-8")
        text_document = strip_accents(text_document.lower())
        return re.findall(self.token_pattern, text_document)


class HashingVectorizer(object):
    """Compute term frequencies vectors using hashed term space

    See the Hashing-trick related papers referenced by John Langford on this
    page to get a grasp on the usefulness of this representation:

      http://hunch.net/~jl/projects/hash_reps/index.html

    dim is the number of buckets, higher dim means lower collision rate but
    also higher memory requirements and higher processing times on the
    resulting tfidf vectors.

    Documents is a sequence of lists of tokens to initialize the DF estimates.

    TODO handle bigrams in a smart way such as demonstrated here:
      http://streamhacker.com/2010/05/24/text-classification-sentiment-analysis-stopwords-collocations/

    """
    # TODO: implement me using the murmurhash that might be faster: but profile
    # me first :)

    def __init__(self, dim=5000, probes=3, analyzer=SimpleAnalyzer()):
        self.dim = dim
        self.probes = probes
        self.analyzer = analyzer

        # start counts at one to avoid zero division while
        # computing IDF
        self.df_counts = np.ones(dim, dtype=long)
        self.tf_vectors = None
        self.labels = []
        self.sampled = 0

    def hash_sign(self, token, probe=0):
        h = hash(token + (probe * u"#"))
        return abs(h) % self.dim, 1.0 if h % 2 == 0 else -1.0

    def sample_document(self, text, tf_vector=None, label=None):
        """Extract features from text and update running freq estimates"""
        if tf_vector is None:
            # allocate term frequency vector and stack to history
            tf_vector = np.zeros(self.dim, np.float64)
            if self.tf_vectors is None:
                self.tf_vectors = tf_vector.reshape((1, self.dim))
            else:
                self.tf_vectors = np.vstack((self.tf_vectors, tf_vector))
                tf_vector = self.tf_vectors[-1]
        if label is not None:
            self.labels.append(label)

        tokens = self.analyzer.analyze(text)
        for token in tokens:
            # TODO add support for cooccurence tokens in a sentence
            # window
            for probe in xrange(self.probes):
                i, incr = self.hash_sign(token, probe)
                tf_vector[i] += incr
        tf_vector /= len(tokens) * self.probes

        # update the running DF estimate
        self.df_counts += tf_vector != 0.0
        self.sampled += 1
        return tf_vector

    def get_tfidf(self):
        """Compute the TF-log(IDF) vectors of the sampled documents"""
        return self.tf_vectors * np.log(float(self.sampled) / self.df_counts)


if __name__ == "__main__":

    # TODO: write unittests instead!
    hv = HashingVectorizer(dim=10, probes=2)
    print hv.analyzer.analyze(u"This is a s\xe9ntence named Mary; with puncts...")

    print hv.sample_document("This is a test document.")
    print hv.sample_document("This is not a test.")
    print hv.sample_document("document document document this toto")
    print hv.tf_vectors
    print hv.df_counts
    print hv.sampled
    print hv.get_tfidf()

