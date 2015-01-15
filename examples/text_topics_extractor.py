"""
=======================================
Extract topics from text using LDA/LSI
=======================================

This is an example showing how scikit-learn can be used extract latent
topic features from a collection of documents. The two topic models used are
Latent Dirichlet Allocation (LDA) and Latent Semantic Indexing (LSI).
The learned topic weights of the documents are saved in CSR sparse matrix
format, and are used as new features representing the original documents.

The dataset used in this example is the 20 newsgroups dataset. It will be
automatically downloaded, then cached.

This example is written by Nan Li (nanli@odesk.com,
nanli@alumni.cs.ucsb.edu).
"""
import sys
sys.path.insert(0, '..')
import logging
from argparse import ArgumentParser

from sklearn.feature_extraction.text import LdaVectorizer, LsiVectorizer
from sklearn.datasets import fetch_20newsgroups

logging.basicConfig(level=logging.INFO)

# fetch the data
newsgroups_train = fetch_20newsgroups(subset='train')

# parse commandline arguments
parser = ArgumentParser()
parser.add_argument('-n', '--num_topics', type=int, default=20)
parser.add_argument('-min_df', '--min_doc_freq', type=float, default=0.001)
parser.add_argument('-max_df', '--max_doc_freq', type=float, default=0.5)
parser.add_argument('-ngram_min', '--ngram_min', type=int, default=1)
parser.add_argument('-ngram_max', '--ngram_max', type=int, default=1)
parser.add_argument('-write_topics', '--write_topics', action='store_true')
parser.add_argument('-lda', '--use_lda', action='store_true')
parser.add_argument('-lsi', '--use_lsi', action='store_true')

usage = '''Usage: text_topics_extractor.py -n <num_topics>
           -min_df <min_doc_freq> -max_df <max_doc_freq>
           -ngram_min <ngram_min> -ngram_max <ngram_max>
           -write_topics (-lda|-lsi)'''
try:
    args = parser.parse_args()
except SystemExit:
    print usage
    sys.exit(2)

ngram_range = args.ngram_min, args.ngram_max

# if LDA model is used
if args.use_lda:
    lda_v = LdaVectorizer(num_topics=args.num_topics,
                          min_df=args.min_doc_freq,
                          max_df=args.max_doc_freq,
                          ngram_range=ngram_range,
                          writedown_topics=args.write_topics)
    lda_weights = lda_v.fit_transform(newsgroups_train.data)

# if LSI model is used
if args.use_lsi:
    lsi_v = LsiVectorizer(num_topics=args.num_topics,
                          min_df=args.min_doc_freq,
                          max_df=args.max_doc_freq,
                          ngram_range=ngram_range,
                          writedown_topics=args.write_topics)
    lsi_weights = lsi_v.fit_transform(newsgroups_train.data)
