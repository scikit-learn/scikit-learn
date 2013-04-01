# -*- coding: utf-8 -*-
# Authors: Olivier Grisel <olivier.grisel@ensta.org>
#          Mathieu Blondel <mathieu@mblondel.org>
#          Lars Buitinck <L.J.Buitinck@uva.nl>
#          Robert Layton <robertlayton@gmail.com>
#          Jochen Wersdörfer <jochen@wersdoerfer.de>
#          Roman Sinayev <roman.sinayev@gmail.com>
#
# License: BSD Style.
"""
The :mod:`sklearn.feature_extraction.text` submodule gathers utilities to
build feature vectors from text documents.
"""
from __future__ import unicode_literals

import array
from collections import Mapping, defaultdict
import numbers
from operator import itemgetter
import re
import unicodedata
import warnings
from multiprocessing import Queue, Process, cpu_count

import numpy as np
import scipy.sparse as sp

from ..base import BaseEstimator, TransformerMixin
from ..externals.six.moves import xrange
from ..preprocessing import normalize
from .hashing import FeatureHasher
from .stop_words import ENGLISH_STOP_WORDS
from sklearn.externals import six

__all__ = ['CountVectorizer',
           'ENGLISH_STOP_WORDS',
           'TfidfTransformer',
           'TfidfVectorizer',
           'strip_accents_ascii',
           'strip_accents_unicode',
           'strip_tags']


def strip_accents_unicode(s):
    """Transform accentuated unicode symbols into their simple counterpart

    Warning: the python-level loop and join operations make this
    implementation 20 times slower than the strip_accents_ascii basic
    normalization.

    See also
    --------
    strip_accents_ascii
        Remove accentuated char for any unicode symbol that has a direct
        ASCII equivalent.
    """
    return ''.join([c for c in unicodedata.normalize('NFKD', s)
                    if not unicodedata.combining(c)])


def strip_accents_ascii(s):
    """Transform accentuated unicode symbols into ascii or nothing

    Warning: this solution is only suited for languages that have a direct
    transliteration to ASCII symbols.

    See also
    --------
    strip_accents_unicode
        Remove accentuated char for any unicode symbol.
    """
    nkfd_form = unicodedata.normalize('NFKD', s)
    return nkfd_form.encode('ASCII', 'ignore').decode('ASCII')


def strip_tags(s):
    """Basic regexp based HTML / XML tag stripper function

    For serious HTML/XML preprocessing you should rather use an external
    library such as lxml or BeautifulSoup.
    """
    return re.compile(r"<([^>]+)>", flags=re.UNICODE).sub(" ", s)


def _check_stop_list(stop):
    if stop == "english":
        return ENGLISH_STOP_WORDS
    elif isinstance(stop, six.string_types):
        raise ValueError("not a built-in stop list: %s" % stop)
    else:               # assume it's a collection
        return stop


class VectorizerMixin(object):
    """Provides common code for text vectorizers (tokenization logic)."""

    _white_spaces = re.compile(r"\s\s+")

    def decode(self, doc):
        """Decode the input into a string of unicode symbols

        The decoding strategy depends on the vectorizer parameters.
        """
        if self.input == 'filename':
            with open(doc, 'rb') as fh:
                doc = fh.read()

        elif self.input == 'file':
            doc = doc.read()

        if isinstance(doc, bytes):
            doc = doc.decode(self.charset, self.charset_error)
        return doc

    def _word_ngrams(self, tokens, stop_words=None):
        """Turn tokens into a sequence of n-grams after stop words filtering"""
        # handle stop words
        if stop_words is not None:
            tokens = [w for w in tokens if w not in stop_words]

        # handle token n-grams
        min_n, max_n = self.ngram_range
        if max_n != 1:
            original_tokens = tokens
            tokens = []
            n_original_tokens = len(original_tokens)
            for n in xrange(min_n,
                            min(max_n + 1, n_original_tokens + 1)):
                for i in xrange(n_original_tokens - n + 1):
                    tokens.append(" ".join(original_tokens[i: i + n]))

        return tokens

    def _char_ngrams(self, text_document):
        """Tokenize text_document into a sequence of character n-grams"""
        # normalize white spaces
        text_document = self._white_spaces.sub(" ", text_document)

        text_len = len(text_document)
        ngrams = []
        min_n, max_n = self.ngram_range
        for n in xrange(min_n, min(max_n + 1, text_len + 1)):
            for i in xrange(text_len - n + 1):
                ngrams.append(text_document[i: i + n])
        return ngrams

    def _char_wb_ngrams(self, text_document):
        """Whitespace sensitive char-n-gram tokenization.

        Tokenize text_document into a sequence of character n-grams
        excluding any whitespace (operating only inside word boundaries)"""
        # normalize white spaces
        text_document = self._white_spaces.sub(u" ", text_document)

        min_n, max_n = self.ngram_range
        ngrams = []
        for w in text_document.split():
            w = ' ' + w + ' '
            w_len = len(w)
            for n in xrange(min_n, max_n + 1):
                offset = 0
                ngrams.append(w[offset:offset + n])
                while offset + n < w_len:
                    offset += 1
                    ngrams.append(w[offset:offset + n])
                if offset == 0:   # count a short word (w_len < n) only once
                    break
        return ngrams

    def build_preprocessor(self):
        """Return a function to preprocess the text before tokenization"""
        if self.preprocessor is not None:
            return self.preprocessor

        # unfortunately python functools package does not have an efficient
        # `compose` function that would have allowed us to chain a dynamic
        # number of functions. However the cost of a lambda call is a few
        # hundreds of nanoseconds which is negligible when compared to the
        # cost of tokenizing a string of 1000 chars for instance.
        noop = lambda x: x

        # accent stripping
        if not self.strip_accents:
            strip_accents = noop
        elif callable(self.strip_accents):
            strip_accents = self.strip_accents
        elif self.strip_accents == 'ascii':
            strip_accents = strip_accents_ascii
        elif self.strip_accents == 'unicode':
            strip_accents = strip_accents_unicode
        else:
            raise ValueError('Invalid value for "strip_accents": %s' %
                             self.strip_accents)

        if self.lowercase:
            return lambda x: strip_accents(x.lower())
        else:
            return strip_accents

    def build_tokenizer(self):
        """Return a function that split a string in sequence of tokens"""
        if self.tokenizer is not None:
            return self.tokenizer
        token_pattern = re.compile(self.token_pattern)
        return lambda doc: token_pattern.findall(doc)

    def get_stop_words(self):
        """Build or fetch the effective stop words list"""
        return _check_stop_list(self.stop_words)

    def build_analyzer(self):
        """Return a callable that handles preprocessing and tokenization"""
        if callable(self.analyzer):
            return self.analyzer

        preprocess = self.build_preprocessor()

        if self.analyzer == 'char':
            return lambda doc: self._char_ngrams(preprocess(self.decode(doc)))

        elif self.analyzer == 'char_wb':
            return lambda doc: self._char_wb_ngrams(
                preprocess(self.decode(doc)))

        elif self.analyzer == 'word':
            stop_words = self.get_stop_words()
            tokenize = self.build_tokenizer()

            return lambda doc: self._word_ngrams(
                tokenize(preprocess(self.decode(doc))), stop_words)

        else:
            raise ValueError('%s is not a valid tokenization scheme/analyzer' %
                             self.analyzer)


class HashingVectorizer(BaseEstimator, VectorizerMixin):
    """Convert a collection of text documents to a matrix of token occurrences

    It turns a collection of text documents into a scipy.sparse matrix holding
    token occurrence counts (or binary occurrence information), possibly
    normalized as token frequencies if norm='l1' or projected on the euclidean
    unit sphere if norm='l2'.

    This text vectorizer implementation uses the hashing trick to find the
    token string name to feature integer index mapping.

    This strategy has several advantage:

    - it is very low memory scalable to large datasets as there is no need to
      store a vocabulary dictionary in memory

    - it is fast to pickle and un-pickle has it holds no state besides the
      constructor parameters

    - it can be used in a streaming (partial fit) or parallel pipeline as there
      is no state computed during fit.

    There are also a couple of cons (vs using a CountVectorizer with an
    in-memory vocabulary):

    - there is no way to compute the inverse transform (from feature indices to
      string feature names) which can be a problem when trying to introspect
      which features are most important to a model.

    - there can be collisions: distinct tokens can be mapped to the same
      feature index. However in practice this is rarely an issue if n_features
      is large enough (e.g. 2 ** 18 for text classification problems).

    - no IDF weighting as this would render the transformer stateful.

    The hash function employed is the signed 32-bit version of Murmurhash3.

    Parameters
    ----------

    input: string {'filename', 'file', 'content'}
        If filename, the sequence passed as an argument to fit is
        expected to be a list of filenames that need reading to fetch
        the raw content to analyze.

        If 'file', the sequence items must have 'read' method (file-like
        object) it is called to fetch the bytes in memory.

        Otherwise the input is expected to be the sequence strings or
        bytes items are expected to be analyzed directly.

    charset: string, 'utf-8' by default.
        If bytes or files are given to analyze, this charset is used to
        decode.

    charset_error: {'strict', 'ignore', 'replace'}
        Instruction on what to do if a byte sequence is given to analyze that
        contains characters not of the given `charset`. By default, it is
        'strict', meaning that a UnicodeDecodeError will be raised. Other
        values are 'ignore' and 'replace'.

    strip_accents: {'ascii', 'unicode', None}
        Remove accents during the preprocessing step.
        'ascii' is a fast method that only works on characters that have
        an direct ASCII mapping.
        'unicode' is a slightly slower method that works on any characters.
        None (default) does nothing.

    analyzer: string, {'word', 'char', 'char_wb'} or callable
        Whether the feature should be made of word or character n-grams.
        Option 'char_wb' creates character n-grams only from text inside
        word boundaries.

        If a callable is passed it is used to extract the sequence of features
        out of the raw, unprocessed input.

    preprocessor: callable or None (default)
        Override the preprocessing (string transformation) stage while
        preserving the tokenizing and n-grams generation steps.

    tokenizer: callable or None (default)
        Override the string tokenization step while preserving the
        preprocessing and n-grams generation steps.

    ngram_range: tuple (min_n, max_n)
        The lower and upper boundary of the range of n-values for different
        n-grams to be extracted. All values of n such that min_n <= n <= max_n
        will be used.

    stop_words: string {'english'}, list, or None (default)
        If a string, it is passed to _check_stop_list and the appropriate stop
        list is returned. 'english' is currently the only supported string
        value.

        If a list, that list is assumed to contain stop words, all of which
        will be removed from the resulting tokens.

    lowercase: boolean, default True
        Convert all characters to lowercase before tokenizing.

    token_pattern: string
        Regular expression denoting what constitutes a "token", only used
        if `tokenize == 'word'`. The default regexp select tokens of 2
        or more letters characters (punctuation is completely ignored
        and always treated as a token separator).

    n_features : interger, optional, (2 ** 20) by default
        The number of features (columns) in the output matrices. Small numbers
        of features are likely to cause hash collisions, but large numbers
        will cause larger coefficient dimensions in linear learners.

    norm : 'l1', 'l2' or None, optional
        Norm used to normalize term vectors. None for no normalization.

    binary: boolean, False by default.
        If True, all non zero counts are set to 1. This is useful for discrete
        probabilistic models that model binary events rather than integer
        counts.

    dtype: type, optional
        Type of the matrix returned by fit_transform() or transform().

    non_negative : boolean, optional
        Whether output matrices should contain non-negative values only;
        effectively calls abs on the matrix prior to returning it.
        When True, output values will be multinomially distributed.
        When False, output values will be normally distributed (Gaussian) with
        mean 0, assuming a good hash function.

    See also
    --------
    CountVectorizer, TfidfVectorizer

    """
    def __init__(self, input='content', charset='utf-8',
                 charset_error='strict', strip_accents=None,
                 lowercase=True, preprocessor=None, tokenizer=None,
                 stop_words=None, token_pattern=r"(?u)\b\w\w+\b",
                 ngram_range=(1, 1), analyzer='word', n_features=(2 ** 20),
                 binary=False, norm='l2', non_negative=False,
                 dtype=np.float64):
        self.input = input
        self.charset = charset
        self.charset_error = charset_error
        self.strip_accents = strip_accents
        self.preprocessor = preprocessor
        self.tokenizer = tokenizer
        self.analyzer = analyzer
        self.lowercase = lowercase
        self.token_pattern = token_pattern
        self.stop_words = stop_words
        self.n_features = n_features
        self.ngram_range = ngram_range
        self.binary = binary
        self.norm = norm
        self.non_negative = non_negative
        self.dtype = dtype

    def partial_fit(self, X, y=None):
        """Does nothing: this transformer is stateless.

        This method is just there to mark the fact that this transformer
        can work in a streaming setup.

        """
        return self

    def fit(self, X, y=None):
        """Does nothing: this transformer is stateless."""
        # triggers a parameter validation
        self._get_hasher().fit(X, y=y)
        return self

    def transform(self, X, y=None):
        """Transform a sequence of instances to a scipy.sparse matrix.

        Parameters
        ----------
        X : iterable over raw text documents, length = n_samples
            Samples. Each sample must be a text document (either bytes or
            unicode strings, filen ame or file object depending on the
            constructor argument) which will be tokenized and hashed.

        y : (ignored)

        Returns
        -------
        X : scipy.sparse matrix, shape = (n_samples, self.n_features)
            Feature matrix, for use with estimators or further transformers.

        """
        analyzer = self.build_analyzer()
        X = self._get_hasher().transform(analyzer(doc) for doc in X)
        if self.binary:
            X.data.fill(1)
        if self.norm is not None:
            X = normalize(X, norm=self.norm, copy=False)
        return X

    # Alias transform to fit_transform for convenience
    fit_transform = transform

    def _get_hasher(self):
        return FeatureHasher(n_features=self.n_features,
                             input_type='string', dtype=self.dtype,
                             non_negative=self.non_negative)


class CountVectorizer(BaseEstimator, VectorizerMixin):
    """Convert a collection of text documents to a matrix of token counts

    This implementation produces a sparse representation of the counts using
    scipy.sparse.coo_matrix.

    If you do not provide an a-priori dictionary and you do not use an analyzer
    that does some kind of feature selection then the number of features will
    be equal to the vocabulary size found by analyzing the data.

    Parameters
    ----------
    input : string {'filename', 'file', 'content'}
        If filename, the sequence passed as an argument to fit is
        expected to be a list of filenames that need reading to fetch
        the raw content to analyze.

        If 'file', the sequence items must have 'read' method (file-like
        object) it is called to fetch the bytes in memory.

        Otherwise the input is expected to be the sequence strings or
        bytes items are expected to be analyzed directly.

    charset : string, 'utf-8' by default.
        If bytes or files are given to analyze, this charset is used to
        decode.

    charset_error : {'strict', 'ignore', 'replace'}
        Instruction on what to do if a byte sequence is given to analyze that
        contains characters not of the given `charset`. By default, it is
        'strict', meaning that a UnicodeDecodeError will be raised. Other
        values are 'ignore' and 'replace'.

    strip_accents : {'ascii', 'unicode', None}
        Remove accents during the preprocessing step.
        'ascii' is a fast method that only works on characters that have
        an direct ASCII mapping.
        'unicode' is a slightly slower method that works on any characters.
        None (default) does nothing.

    analyzer : string, {'word', 'char', 'char_wb'} or callable
        Whether the feature should be made of word or character n-grams.
        Option 'char_wb' creates character n-grams only from text inside
        word boundaries.

        If a callable is passed it is used to extract the sequence of features
        out of the raw, unprocessed input.

    preprocessor : callable or None (default)
        Override the preprocessing (string transformation) stage while
        preserving the tokenizing and n-grams generation steps.

    tokenizer : callable or None (default)
        Override the string tokenization step while preserving the
        preprocessing and n-grams generation steps.

    ngram_range : tuple (min_n, max_n)
        The lower and upper boundary of the range of n-values for different
        n-grams to be extracted. All values of n such that min_n <= n <= max_n
        will be used.

    stop_words : string {'english'}, list, or None (default)
        If a string, it is passed to _check_stop_list and the appropriate stop
        list is returned. 'english' is currently the only supported string
        value.

        If a list, that list is assumed to contain stop words, all of which
        will be removed from the resulting tokens.

        If None, no stop words will be used. max_df can be set to a value
        in the range [0.7, 1.0) to automatically detect and filter stop
        words based on intra corpus document frequency of terms.

    lowercase : boolean, default True
        Convert all characters to lowercase befor tokenizing.

    token_pattern : string
        Regular expression denoting what constitutes a "token", only used
        if `tokenize == 'word'`. The default regexp select tokens of 2
        or more letters characters (punctuation is completely ignored
        and always treated as a token separator).

    max_df : float in range [0.0, 1.0] or int, optional, 1.0 by default
        When building the vocabulary ignore terms that have a term frequency
        strictly higher than the given threshold (corpus specific stop words).
        If float, the parameter represents a proportion of documents, integer
        absolute counts.
        This parameter is ignored if vocabulary is not None.

    min_df : float in range [0.0, 1.0] or int, optional, 1 by default
        When building the vocabulary ignore terms that have a term frequency
        strictly lower than the given threshold. This value is also called
        cut-off in the literature.
        If float, the parameter represents a proportion of documents, integer
        absolute counts.
        This parameter is ignored if vocabulary is not None.

    max_features : optional, None by default
        If not None, build a vocabulary that only consider the top
        max_features ordered by term frequency across the corpus.

        This parameter is ignored if vocabulary is not None.

    vocabulary : Mapping or iterable, optional
        Either a Mapping (e.g., a dict) where keys are terms and values are
        indices in the feature matrix, or an iterable over terms. If not
        given, a vocabulary is determined from the input documents.

    binary : boolean, False by default.
        If True, all non zero counts are set to 1. This is useful for discrete
        probabilistic models that model binary events rather than integer
        counts.

    dtype : type, optional
        Type of the matrix returned by fit_transform() or transform().

    sort_features : boolean, True by default.
        If True, the output matrix has the feature columns in order
        corresponding to sorted vocabulary.  If False, order of the
        columns is determined by features first encountered.
        Always True if n_jobs > 1.


    Attributes
    ----------
    `vocabulary_` : dict
        A mapping of terms to feature indices.

    `stop_words_` : set
        Terms that were ignored because
        they occurred in either too many
        (`max_df`) or in too few (`min_df`) documents.
        This is only available if no vocabulary was given.

    See also
    --------
    HashingVectorizer, TfidfVectorizer
    """

    def __init__(self, input='content', charset='utf-8',
                 charset_error='strict', strip_accents=None,
                 lowercase=True, preprocessor=None, tokenizer=None,
                 stop_words=None, token_pattern=r"(?u)\b\w\w+\b",
                 ngram_range=(1, 1), analyzer='word',
                 max_df=1.0, min_df=1, max_features=None,
                 vocabulary=None, binary=False, dtype=np.int64,
                 sort_features=True, n_jobs=1):
        self.input = input
        self.charset = charset
        self.charset_error = charset_error
        self.strip_accents = strip_accents
        self.preprocessor = preprocessor
        self.tokenizer = tokenizer
        self.analyzer = analyzer
        self.lowercase = lowercase
        self.token_pattern = token_pattern
        self.stop_words = stop_words
        self.max_df = max_df
        self.min_df = min_df
        self.n_jobs = n_jobs
        if max_df < 0 or min_df < 0:
            raise ValueError("negative value for max_df of min_df")
        self.max_features = max_features
        if not any((
                isinstance(max_features, numbers.Integral),
                max_features is None,
                max_features > 0)):
            raise ValueError(
                "max_features is neither a positive integer nor None")
        self.ngram_range = ngram_range
        if vocabulary is not None:
            if not isinstance(vocabulary, Mapping):
                vocabulary = dict((t, i) for i, t in enumerate(vocabulary))
            if not vocabulary:
                raise ValueError("empty vocabulary passed to fit")
            self.fixed_vocabulary = True
            self.vocabulary_ = dict(vocabulary)
        else:
            self.fixed_vocabulary = False
        self.binary = binary
        self.dtype = dtype
        self.sort_features = sort_features

    def _term_counts_to_matrix(self, n_doc, i_indices,
                               j_indices, values, n_features):
        """Construct COO matrix from indices and values."""
        # j_indices is either an array.array or numpy array of np.int32 dtype
        if type(j_indices) is array.array:
            try:
                j_indices = np.frombuffer(j_indices, dtype=np.intc)
            except ValueError:
                j_indices = np.array([])
        shape = (n_doc, n_features)
        spmatrix = sp.coo_matrix((values, (i_indices, j_indices)),
                                 shape=shape, dtype=self.dtype)
        # weird bug -- this doesn't work on some computers
        # if self.binary:
        #     spmatrix.data.fill(1)
        #     print "self.binary", self.binary
        #     print spmatrix.todense()
        return spmatrix

    def _sort_features_and_matrix(self, cscmatrix, feature_to_position):
        '''sort dict by keys and assign values to the sorted key index'''
        sorted_by_name_dict = {}
        sorted_features = sorted(feature_to_position)
        new_positions = []
        for i, k in enumerate(sorted_features):
            sorted_by_name_dict[k] = i
            new_positions.append(feature_to_position[k])
        return cscmatrix[:, new_positions], sorted_by_name_dict

    def _remove_highandlow(self, cscmatrix,
                           feature_to_position, high, low):
        """Remove too rare or too common features.

        Prune features that are non zero in more samples than high or less
        documents than low.

        This does not prune samples with zero features.

        """
        kept_indices = []
        removed_indices = set()
        for colptr in xrange(len(cscmatrix.indptr) - 1):
            len_data_slice = len(cscmatrix.data[cscmatrix.indptr[colptr]:
                                                cscmatrix.indptr[colptr + 1]])
            if len_data_slice <= high and len_data_slice >= low:
                kept_indices.append(colptr)
            else:
                removed_indices.add(colptr)
        s_kept_indices = set(kept_indices)
        new_mapping = dict((v, i) for i, v in enumerate(kept_indices))
        feature_to_position = dict((k, new_mapping[v]) for k, v in
                                   six.iteritems(feature_to_position)
                                   if v in s_kept_indices)
        return cscmatrix[:, kept_indices], feature_to_position, removed_indices

    def _get_kept_features(self, csc_m, max_features):
        '''Helper method for _get_max_features to use less memory.'''
        feature_freqs = np.asarray(csc_m.sum(axis=0)).ravel()
        to_keep = np.argsort(feature_freqs)[::-1][:max_features]
        return to_keep

    def _get_max_features(self, csc_m, feature_to_position,
                          stop_words_, max_features):
        '''Remove maximum features using a sparse matrix.

        Cut only the top max_features from the matrix and
        feature_to_position. Cut features are added to stop_words_.

        '''
        to_keep = self._get_kept_features(csc_m, max_features)
        s_to_keep = set(to_keep)
        new_feature_to_position = {}
        for (k, v) in six.iteritems(feature_to_position):
            if v in s_to_keep:
                new_feature_to_position[k] = v
            else:
                stop_words_.add(k)
        return csc_m[:, np.sort(to_keep)], new_feature_to_position, stop_words_

    def _merge_outputs(self, to_merge_list):
        '''Helper method for _parallel_count. Merge _count_*_vocab outputs.'''
        j_ind_old, f2p_old = to_merge_list[0]
        for (j_ind, f2p) in to_merge_list[1:]:
            # remap j_indices to account for words previously encountered
            # by other processes
            position_remap = {}
            last_feature_index = len(f2p_old) - 1
            for key in f2p:
                if key in f2p_old:
                    position_remap[f2p[key]] = f2p_old[key]
                else:
                    last_feature_index += 1
                    # update the feature dict
                    f2p_old[key] = last_feature_index
                    position_remap[f2p[key]] = last_feature_index
            # concatenate j_ind and j_ind_old according to the remapping.
            # the below loop takes 1/3 of the time of the whole fit_transform.
            # needs to be optimized.
            # both j_ind and j_ind_old are "i" array.array objects
            j_ind_old.extend([position_remap[j] for j in j_ind])
        return j_ind_old, f2p_old

    def _run_queue(self, fixed_vocab, in_q, out_q):
        """Run the queue created in _parallel_count."""
        idx, chunk = in_q.get()
        if fixed_vocab:
            out_q.put((idx, self._count_fixed_vocab(chunk)))
        else:
            out_q.put((idx, self._count_new_vocab(chunk)))

    def _parallel_count(self, raw_documents, fixed_vocab, n_jobs):
        """Use mp to parallelize _count_new_vocab and _count_fixed_vocab.

        Chunk raw documents into n_jobs lists and apply _count_*_vocab on each
        using pythone Queue module.
        Then merge the outputs using _merge_outputs.

        """
        if not hasattr(raw_documents, "__getitem__"):
            raw_documents = [doc for doc in raw_documents]
        if len(raw_documents) == 0:
            raise ValueError("No documents provided")
        chunk_size = int(np.ceil(float(len(raw_documents))/n_jobs))
        in_q = Queue(n_jobs)
        out_q = Queue(n_jobs)
        procs = []
        for i, chunk in enumerate(self._chunk_docs(raw_documents, chunk_size)):
            in_q.put((i, chunk))
        n_chunks = i + 1
        for i in xrange(n_chunks):
            procs.append(Process(target=self._run_queue,
                                 args=(fixed_vocab, in_q, out_q)))
        for proc in procs:
            proc.start()
        process_output = [out_q.get() for i in xrange(len(procs))]
        for proc in procs:
            proc.join()
        process_output.sort()
        i_arrays = []
        to_merge_list = []
        total_n_doc = 0
        if not fixed_vocab:
            for i, (j_indices, feature_to_position,
                    n_doc, features_per_doc) in process_output:
                i_array = self._make_i_indices(features_per_doc, total_n_doc)
                total_n_doc += n_doc
                # if i_array is None, no features were found
                # in this process's chunk
                if i_array is None:
                    continue
                i_arrays.append(i_array)
                to_merge_list.append((j_indices, feature_to_position))

            if len(i_arrays) == 0:
                error_str = ''.join(["empty vocabulary;",
                                     " Perhaps the document",
                                     " only contains stop words"])
                raise ValueError(error_str)
            j_indices, feature_to_position = self._merge_outputs(to_merge_list)
            i_indices = np.concatenate(i_arrays)
            j_indices = np.frombuffer(j_indices, dtype=np.intc)
            return i_indices, j_indices, feature_to_position, total_n_doc
        else:  # fixed vocabulary
            for i, (j_indices, n_doc, features_per_doc) in process_output:
                i_array = self._make_i_indices(features_per_doc, total_n_doc)
                total_n_doc += n_doc
                if i_array is None:
                    continue
                i_arrays.append(i_array)
                to_merge_list.append(np.frombuffer(j_indices, dtype=np.intc))
            if len(i_arrays) == 0:
                error_str = ''.join(["empty vocabulary;",
                                     " Perhaps the document",
                                     " only contains stop words"])
                raise ValueError(error_str)
            j_indices = np.concatenate(to_merge_list)
            i_indices = np.concatenate(i_arrays)
            return i_indices, j_indices, total_n_doc

    def _chunk_docs(self, iterable, chunk_size):
        """Chunk iterable docs into chunks of chunk_size."""
        iterable = iter(iterable)
        while True:
            chunk = []
            try:
                for _ in range(chunk_size):
                    chunk.append(iterable.next())
                yield chunk
            except StopIteration:
                if chunk:
                    yield chunk
                break

    def _count_fixed_vocab(self, raw_documents):
        """Make feature position indices, count number of features per doc.

        Follow the same strategy as _count_new_vocab but with known
        feature_to_position.

        """
        j_indices = _make_int_array()
        analyze = self.build_analyzer()
        features_per_doc = []
        for i, doc in enumerate(raw_documents):
            k = 0
            for feature in analyze(doc):
                try:
                    j_indices.append(self.vocabulary_[feature])
                    k += 1
                except KeyError:
                    continue
            features_per_doc.append(k)
        n_doc = i + 1
        return j_indices, n_doc, features_per_doc

    def _count_new_vocab(self, raw_documents):
        """Create feature position indices and idx to matrix column mapping.

        To create a COO matrix, we need three arrays: i_indices, j_indices and
        values.  When a new feature is encountered, we add a number to each.
        Values can all be 1, and i indices are all the same for all features
        within one document, so both can be quickly created later.
        We take advantage of the fact that values in
        the same position get implicitly added together when the COO
        matrix is constructed when i indices and j indices are the same.
        Here we create j_indices and count how many i_indices we need per doc.
        feature_to_position is a mapping between features as strings
        and the final location in the matrix.

        """
        analyze = self.build_analyzer()
        j = 0  # counts new features
        feature_to_count = defaultdict(int)
        feature_to_position = {}
        j_indices = _make_int_array()
        features_per_doc = []
        for i, doc in enumerate(raw_documents):
            for k, feature in enumerate(analyze(doc)):
                feature_to_count[feature] += 1
                if feature_to_count[feature] == 1:  # new feature
                    feature_to_position[feature] = j
                    j += 1
                j_indices.append(feature_to_position[feature])
            try:
                features_per_doc.append(k+1)
            except UnboundLocalError:
                # k is undefined so
                # no features found in the doc
                features_per_doc.append(0)
        # assert j == len(feature_to_count)
        n_doc = i + 1
        return j_indices, feature_to_position, n_doc, features_per_doc

    def _make_i_indices(self, features_per_doc, start_index=0):
        '''Create i_indices from features_per_doc.

        features_per_doc is a list of the number of kept features in each doc
        in order.  When using multiprocessing, each process's docs index may
        not start at zero, so we start counting with start_index instead.

        '''
        combined_arrays = []
        number_of_docs_with_features = 0
        for i, num_features in enumerate(features_per_doc):
            if num_features > 0:  # some features found in the doc
                number_of_docs_with_features += 1
                i_vals = np.empty(num_features, dtype=np.int32)
                fill_val = i + start_index
                i_vals.fill(fill_val)
                combined_arrays.append(i_vals)
        if number_of_docs_with_features == 1:
            return combined_arrays[0]
        elif number_of_docs_with_features == 0:  # no docs found with features
            return None
        else:  # more than one doc with features
            return np.concatenate(combined_arrays)

    def fit(self, raw_documents, y=None):
        """Learn a vocabulary dictionary of all tokens in the raw documents.

        Parameters
        ----------
        raw_documents : iterable
            An iterable which yields either str, unicode or file objects.

        Returns
        -------
        self
        """
        self.fit_transform(raw_documents)
        return self

    def fit_transform(self, raw_documents, y=None):
        """Learn the vocabulary dictionary and return the count vectors.

        This is more efficient than calling fit followed by transform.

        Parameters
        ----------
        raw_documents : iterable
            An iterable which yields either str, unicode or file objects.

        Returns
        -------
        vectors : array, [n_samples, n_features]
        """
        # We intentionally don't call the transform method to make
        # fit_transform overridable without unwanted side effects in
        # TfidfVectorizer.
        fixed_vocab = self.fixed_vocabulary
        max_df = self.max_df
        min_df = self.min_df
        max_features = self.max_features
        sort_features = self.sort_features
        binary = self.binary
        n_jobs = self.n_jobs
        if n_jobs == -1:
            n_jobs = cpu_count()
        if n_jobs == 1:  # use serial execution
            if fixed_vocab:
                feature_to_position = self.vocabulary_
                j_indices, n_doc, features_per_doc = \
                    self._count_fixed_vocab(raw_documents)
                i_indices = self._make_i_indices(features_per_doc)
            else:
                j_indices, feature_to_position, n_doc, features_per_doc = \
                    self._count_new_vocab(raw_documents)
                i_indices = self._make_i_indices(features_per_doc)
            del features_per_doc  # free memory
        else:  # use parallel execution
            if fixed_vocab:
                feature_to_position = self.vocabulary_
                i_indices, j_indices, n_doc = \
                    self._parallel_count(raw_documents, fixed_vocab, n_jobs)
            else:
                i_indices, j_indices, feature_to_position, n_doc = \
                    self._parallel_count(raw_documents, fixed_vocab, n_jobs)
        if not fixed_vocab:
            max_doc_count = (max_df
                             if isinstance(max_df, numbers.Integral)
                             else int(round(max_df * n_doc))
                             )
            min_doc_count = (min_df
                             if isinstance(min_df, numbers.Integral)
                             else int(round(min_df * n_doc))
                             )
            if max_doc_count < min_doc_count:
                raise ValueError(
                    "max_df corresponds to < documents than min_df")
        n_features = len(feature_to_position)
        if n_features == 0 or i_indices is None:
            error_str = ''.join(["empty vocabulary;",
                                 " Perhaps the document",
                                 " only contains stop words"])
            raise ValueError(error_str)
        values = np.empty(len(i_indices), dtype=np.int32)
        values.fill(1)
        csc_m = self._term_counts_to_matrix(n_doc, i_indices, j_indices,
                                            values, n_features).tocsc()
        del i_indices, j_indices, values  # free memory
        if binary:
            csc_m.data.fill(1)
        if not fixed_vocab:
            if sort_features or n_jobs != 1:
                csc_m, feature_to_position = \
                    self._sort_features_and_matrix(csc_m, feature_to_position)
            stop_words_ = set()
            # get rid of features between max_df and min_df
            if max_doc_count < n_doc or min_doc_count > 1:
                csc_m, feature_to_position, stop_words_ = \
                    self._remove_highandlow(csc_m, feature_to_position,
                                            max_doc_count, min_doc_count
                                            )
            # get rid of features that are not in the top max_features
            # overall occurance wise
            if max_features and max_features < csc_m.shape[1]:
                csc_m, feature_to_position, stop_words_ = \
                    self._get_max_features(
                        csc_m, feature_to_position, stop_words_, max_features)
            self.stop_words_ = stop_words_
            self.vocabulary_ = feature_to_position
        return csc_m

    def transform(self, raw_documents):
        """Extract token counts out of raw text documents using the vocabulary
        fitted with fit or the one provided in the constructor.

        Parameters
        ----------
        raw_documents : iterable
            An iterable which yields either str, unicode or file objects.

        Returns
        -------
        vectors : sparse matrix, [n_samples, n_features]
        """
        if not hasattr(self, 'vocabulary_') or len(self.vocabulary_) == 0:
            raise ValueError("Vocabulary wasn't fitted or is empty!")

        # raw_documents can be an iterable so we don't know its size in
        # advance
        j_indices = _make_int_array()
        binary = self.binary
        n_jobs = self.n_jobs
        if n_jobs == -1:
            n_jobs = cpu_count()
        # use the same matrix-building strategy as fit_transform
        if n_jobs == 1:
            j_indices, n_doc, features_per_doc = \
                self._count_fixed_vocab(raw_documents)
            i_indices = self._make_i_indices(features_per_doc)
        else:
            i_indices, j_indices, n_doc = \
                self._parallel_count(raw_documents, fixed_vocab=True,
                                     n_jobs=n_jobs
                                     )
        if i_indices is None:
            i_indices = np.empty(0, dtype=np.int32)
        values = np.empty(len(i_indices), dtype=np.int32)
        values.fill(1)
        n_features = len(self.vocabulary_)
        m = self._term_counts_to_matrix(n_doc, i_indices,
                                        j_indices, values, n_features)
        del i_indices, j_indices, values  # free memory
        if binary:
            m.data.fill(1)
        return m.tocsc()  # keep the same sparse format as fit_transform

    def inverse_transform(self, X):
        """Return terms per document with nonzero entries in X.

        Parameters
        ----------
        X : {array, sparse matrix}, shape = [n_samples, n_features]

        Returns
        -------
        X_inv : list of arrays, len = n_samples
            List of arrays of terms.
        """
        if sp.isspmatrix_coo(X) or sp.isspmatrix_csc(X):
            # COO matrix is not indexable, CSC is slow for row manipulations
            X = X.tocsr()
        elif not sp.issparse(X):
            # We need to convert X to a matrix, so that the indexing
            # returns 2D objects
            X = np.asmatrix(X)
        n_samples = X.shape[0]

        terms = np.array(list(self.vocabulary_.keys()))
        indices = np.array(list(self.vocabulary_.values()))
        inverse_vocabulary = terms[np.argsort(indices)]

        return [inverse_vocabulary[X[i, :].nonzero()[1]].ravel()
                for i in range(n_samples)]

    def get_feature_names(self):
        """Array mapping from feature integer indices to feature name"""
        if not hasattr(self, 'vocabulary_') or len(self.vocabulary_) == 0:
            raise ValueError("Vocabulary wasn't fitted or is empty!")

        return [t for t, i in sorted(six.iteritems(self.vocabulary_),
                                     key=itemgetter(1))]

    @property
    def max_df_stop_words_(self):
        warnings.warn(
            "The 'stop_words_ attribute was renamed to 'max_df_stop_words'. "
            "The old attribute will be removed in 0.15.", DeprecationWarning)
        return self.stop_words_


def _make_int_array():
    """Construct an array.array of a type suitable for scipy.sparse indices."""
    return array.array(str("i"))


class TfidfTransformer(BaseEstimator, TransformerMixin):
    """Transform a count matrix to a normalized tf or tf–idf representation

    Tf means term-frequency while tf–idf means term-frequency times inverse
    document-frequency. This is a common term weighting scheme in information
    retrieval, that has also found good use in document classification.

    The goal of using tf–idf instead of the raw frequencies of occurrence of a
    token in a given document is to scale down the impact of tokens that occur
    very frequently in a given corpus and that are hence empirically less
    informative than features that occur in a small fraction of the training
    corpus.

    In the SMART notation used in IR, this class implements several tf–idf
    variants:

    Tf is "n" (natural) by default, "l" (logarithmic) when sublinear_tf=True.
    Idf is "t" idf is "t" when use_idf is given, "n" (none) otherwise.
    Normalization is "c" (cosine) when norm='l2', "n" (none) when norm=None.

    Parameters
    ----------
    norm : 'l1', 'l2' or None, optional
        Norm used to normalize term vectors. None for no normalization.

    use_idf : boolean, optional
        Enable inverse-document-frequency reweighting.

    smooth_idf : boolean, optional
        Smooth idf weights by adding one to document frequencies, as if an
        extra document was seen containing every term in the collection
        exactly once. Prevents zero divisions.

    sublinear_tf : boolean, optional
        Apply sublinear tf scaling, i.e. replace tf with 1 + log(tf).

    References
    ----------

    .. [Yates2011] `R. Baeza-Yates and B. Ribeiro-Neto (2011). Modern
                   Information Retrieval. Addison Wesley, pp. 68–74.`

    .. [MSR2008] `C.D. Manning, H. Schütze and P. Raghavan (2008). Introduction
                 to Information Retrieval. Cambridge University Press,
                 pp. 121–125.`
    """

    def __init__(self, norm='l2', use_idf=True, smooth_idf=True,
                 sublinear_tf=False):
        self.norm = norm
        self.use_idf = use_idf
        self.smooth_idf = smooth_idf
        self.sublinear_tf = sublinear_tf

    def fit(self, X, y=None):
        """Learn the idf vector (global term weights)

        Parameters
        ----------
        X : sparse matrix, [n_samples, n_features]
            a matrix of term/token counts
        """
        if self.use_idf:
            if not hasattr(X, 'nonzero'):
                X = sp.csr_matrix(X)

            n_samples, n_features = X.shape
            df = np.bincount(X.nonzero()[1])
            if df.shape[0] < n_features:
                # bincount might return fewer bins than there are features
                df = np.concatenate([df, np.zeros(n_features - df.shape[0])])

            # perform idf smoothing if required
            df += int(self.smooth_idf)
            n_samples += int(self.smooth_idf)

            # avoid division by zeros for features that occur in all documents
            idf = np.log(float(n_samples) / df) + 1.0
            idf_diag = sp.lil_matrix((n_features, n_features))
            idf_diag.setdiag(idf)
            self._idf_diag = sp.csc_matrix(idf_diag)

        return self

    def transform(self, X, copy=True):
        """Transform a count matrix to a tf or tf–idf representation

        Parameters
        ----------
        X : sparse matrix, [n_samples, n_features]
            a matrix of term/token counts

        Returns
        -------
        vectors : sparse matrix, [n_samples, n_features]
        """
        if hasattr(X, 'dtype') and np.issubdtype(X.dtype, np.float):
            # preserve float family dtype
            X = sp.csr_matrix(X, copy=copy)
        else:
            # convert counts or binary occurrences to floats
            X = sp.csr_matrix(X, dtype=np.float64, copy=copy)

        n_samples, n_features = X.shape

        if self.sublinear_tf:
            np.log(X.data, X.data)
            X.data += 1

        if self.use_idf:
            if not hasattr(self, "_idf_diag"):
                raise ValueError("idf vector not fitted")
            expected_n_features = self._idf_diag.shape[0]
            if n_features != expected_n_features:
                raise ValueError("Input has n_features=%d while the model"
                                 " has been trained with n_features=%d" % (
                                     n_features, expected_n_features))
            # *= doesn't work
            X = X * self._idf_diag

        if self.norm:
            X = normalize(X, norm=self.norm, copy=False)

        return X

    @property
    def idf_(self):
        if hasattr(self, "_idf_diag"):
            return np.ravel(self._idf_diag.sum(axis=0))
        else:
            return None


class TfidfVectorizer(CountVectorizer):
    """Convert a collection of raw documents to a matrix of TF-IDF features.

    Equivalent to CountVectorizer followed by TfidfTransformer.

    Parameters
    ----------
    input : string {'filename', 'file', 'content'}
        If filename, the sequence passed as an argument to fit is
        expected to be a list of filenames that need reading to fetch
        the raw content to analyze.

        If 'file', the sequence items must have 'read' method (file-like
        object) it is called to fetch the bytes in memory.

        Otherwise the input is expected to be the sequence strings or
        bytes items are expected to be analyzed directly.

    charset : string, 'utf-8' by default.
        If bytes or files are given to analyze, this charset is used to
        decode.

    charset_error : {'strict', 'ignore', 'replace'}
        Instruction on what to do if a byte sequence is given to analyze that
        contains characters not of the given `charset`. By default, it is
        'strict', meaning that a UnicodeDecodeError will be raised. Other
        values are 'ignore' and 'replace'.

    strip_accents : {'ascii', 'unicode', None}
        Remove accents during the preprocessing step.
        'ascii' is a fast method that only works on characters that have
        an direct ASCII mapping.
        'unicode' is a slightly slower method that works on any characters.
        None (default) does nothing.

    analyzer : string, {'word', 'char'} or callable
        Whether the feature should be made of word or character n-grams.

        If a callable is passed it is used to extract the sequence of features
        out of the raw, unprocessed input.

    preprocessor : callable or None (default)
        Override the preprocessing (string transformation) stage while
        preserving the tokenizing and n-grams generation steps.

    tokenizer : callable or None (default)
        Override the string tokenization step while preserving the
        preprocessing and n-grams generation steps.

    ngram_range : tuple (min_n, max_n)
        The lower and upper boundary of the range of n-values for different
        n-grams to be extracted. All values of n such that min_n <= n <= max_n
        will be used.

    stop_words : string {'english'}, list, or None (default)
        If a string, it is passed to _check_stop_list and the appropriate stop
        list is returned. 'english' is currently the only supported string
        value.

        If a list, that list is assumed to contain stop words, all of which
        will be removed from the resulting tokens.

        If None, no stop words will be used. max_df can be set to a value
        in the range [0.7, 1.0) to automatically detect and filter stop
        words based on intra corpus document frequency of terms.

    lowercase : boolean, default True
        Convert all characters to lowercase befor tokenizing.

    token_pattern : string
        Regular expression denoting what constitutes a "token", only used
        if `tokenize == 'word'`. The default regexp select tokens of 2
        or more letters characters (punctuation is completely ignored
        and always treated as a token separator).

    max_df : float in range [0.0, 1.0] or int, optional, 1.0 by default
        When building the vocabulary ignore terms that have a term frequency
        strictly higher than the given threshold (corpus specific stop words).
        If float, the parameter represents a proportion of documents, integer
        absolute counts.
        This parameter is ignored if vocabulary is not None.

    min_df : float in range [0.0, 1.0] or int, optional, 1 by default
        When building the vocabulary ignore terms that have a term frequency
        strictly lower than the given threshold.
        This value is also called cut-off in the literature.
        If float, the parameter represents a proportion of documents, integer
        absolute counts.
        This parameter is ignored if vocabulary is not None.

    max_features : optional, None by default
        If not None, build a vocabulary that only consider the top
        max_features ordered by term frequency across the corpus.

        This parameter is ignored if vocabulary is not None.

    vocabulary : Mapping or iterable, optional
        Either a Mapping (e.g., a dict) where keys are terms and values are
        indices in the feature matrix, or an iterable over terms. If not
        given, a vocabulary is determined from the input documents.

    binary : boolean, False by default.
        If True, all non zero counts are set to 1. This is useful for discrete
        probabilistic models that model binary events rather than integer
        counts.

    dtype : type, optional
        Type of the matrix returned by fit_transform() or transform().

    norm : 'l1', 'l2' or None, optional
        Norm used to normalize term vectors. None for no normalization.

    use_idf : boolean, optional
        Enable inverse-document-frequency reweighting.

    smooth_idf : boolean, optional
        Smooth idf weights by adding one to document frequencies, as if an
        extra document was seen containing every term in the collection
        exactly once. Prevents zero divisions.

    sublinear_tf : boolean, optional
        Apply sublinear tf scaling, i.e. replace tf with 1 + log(tf).

    See also
    --------
    CountVectorizer
        Tokenize the documents and count the occurrences of token and return
        them as a sparse matrix

    TfidfTransformer
        Apply Term Frequency Inverse Document Frequency normalization to a
        sparse matrix of occurrence counts.

    """

    def __init__(self, input='content', charset='utf-8',
                 charset_error='strict', strip_accents=None, lowercase=True,
                 preprocessor=None, tokenizer=None, analyzer='word',
                 stop_words=None, token_pattern=r"(?u)\b\w\w+\b",
                 ngram_range=(1, 1), max_df=1.0, min_df=1,
                 max_features=None, vocabulary=None, binary=False,
                 dtype=np.int64, norm='l2', use_idf=True, smooth_idf=True,
                 sublinear_tf=False):

        super(TfidfVectorizer, self).__init__(
            input=input, charset=charset, charset_error=charset_error,
            strip_accents=strip_accents, lowercase=lowercase,
            preprocessor=preprocessor, tokenizer=tokenizer, analyzer=analyzer,
            stop_words=stop_words, token_pattern=token_pattern,
            ngram_range=ngram_range, max_df=max_df, min_df=min_df,
            max_features=max_features, vocabulary=vocabulary, binary=False,
            dtype=dtype)

        self._tfidf = TfidfTransformer(norm=norm, use_idf=use_idf,
                                       smooth_idf=smooth_idf,
                                       sublinear_tf=sublinear_tf)

    # Broadcast the TF-IDF parameters to the underlying transformer instance
    # for easy grid search and repr

    @property
    def norm(self):
        return self._tfidf.norm

    @norm.setter
    def norm(self, value):
        self._tfidf.norm = value

    @property
    def use_idf(self):
        return self._tfidf.use_idf

    @use_idf.setter
    def use_idf(self, value):
        self._tfidf.use_idf = value

    @property
    def smooth_idf(self):
        return self._tfidf.smooth_idf

    @smooth_idf.setter
    def smooth_idf(self, value):
        self._tfidf.smooth_idf = value

    @property
    def sublinear_tf(self):
        return self._tfidf.sublinear_tf

    @sublinear_tf.setter
    def sublinear_tf(self, value):
        self._tfidf.sublinear_tf = value

    def fit(self, raw_documents, y=None):
        """Learn a conversion law from documents to array data"""
        X = super(TfidfVectorizer, self).fit_transform(raw_documents)
        self._tfidf.fit(X)
        return self

    def fit_transform(self, raw_documents, y=None):
        """Learn the representation and return the vectors.

        Parameters
        ----------
        raw_documents : iterable
            an iterable which yields either str, unicode or file objects

        Returns
        -------
        vectors : array, [n_samples, n_features]
        """
        X = super(TfidfVectorizer, self).fit_transform(raw_documents)
        self._tfidf.fit(X)
        # X is already a transformed view of raw_documents so
        # we set copy to False
        return self._tfidf.transform(X, copy=False)

    def transform(self, raw_documents, copy=True):
        """Transform raw text documents to tf–idf vectors

        Parameters
        ----------
        raw_documents : iterable
            an iterable which yields either str, unicode or file objects

        Returns
        -------
        vectors : sparse matrix, [n_samples, n_features]
        """
        X = super(TfidfVectorizer, self).transform(raw_documents)
        return self._tfidf.transform(X, copy)
