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
from argparse import ArgumentParser
sys.path.insert(0, '..')

from sklearn.feature_extraction.text import LdaVectorizer, LsiVectorizer
from sklearn.datasets import fetch_20newsgroups

# fetch the data
newsgroups_train = fetch_20newsgroups(subset='train')

# default parameter values
num_topics = 20
min_df = 10
max_df = 0.5
ngram_min = 1
ngram_max = 1
log_file = ''

# parse commandline arguments
parser = ArgumentParser()
parser.add_argument('-n', '--num_topics')
parser.add_argument('-min_df', '--min_doc_freq')
parser.add_argument('-max_df', '--max_doc_freq')
parser.add_argument('-ngram_min', '--ngram_min')
parser.add_argument('-ngram_max', '--ngram_max')
parser.add_argument('-log_file', '--log_file_name')
parser.add_argument('-lda', '--use_lda', action='store_true')
parser.add_argument('-lsi', '--use_lsi', action='store_true')

usage = '''Usage: text_topics_extractor.py -n <num_topics>
           -min_df <min_doc_freq> -max_df <max_doc_freq>
           -ngram_min <ngram_min> -ngram_max <ngram_max>
           -log_file <log_file_name> (-lda|-lsi)'''
try:
    args = parser.parse_args()
except SystemExit:
    print usage
    sys.exit(2)

if args.num_topics is not None:
    num_topics = int(args.num_topics)
if args.min_doc_freq is not None:
    min_df_f = float(args.min_doc_freq)
    min_df_i = int(min_df_f)
    if min_df_f == min_df_i:
        min_df = min_df_i
    else:
        min_df = min_df_f
if args.max_doc_freq is not None:
    max_df_f = float(args.max_doc_freq)
    max_df_i = int(max_df_f)
    if max_df_f == max_df_i:
        max_df = max_df_i
    else:
        max_df = max_df_f
if args.ngram_min is not None:
    ngram_min = int(args.ngram_min)
if args.ngram_max is not None:
    ngram_max = int(args.ngram_max)
if args.log_file_name is not None:
    log_file = args.log_file_name
elif args.use_lda:
    log_file = 'lda_topics.log'
elif args.use_lsi:
    log_file = 'lsi_topics.log'

ngram_range = ngram_min, ngram_max

# if LDA model is used
if args.use_lda:
    lda_v = LdaVectorizer(num_topics=num_topics,
                          min_df=min_df, max_df=max_df,
                          ngram_range=ngram_range,
                          log_file=log_file)
    lda_weights = lda_v.fit_transform(newsgroups_train.data)

# if LSI model is used
if args.use_lsi:
    lsi_v = LsiVectorizer(num_topics=num_topics,
                          min_df=min_df, max_df=max_df,
                          ngram_range=ngram_range,
                          log_file=log_file)
    lsi_weights = lsi_v.fit_transform(newsgroups_train.data)
