"""
==================================================
Topics extraction with Latent Dirichlet Allocation
==================================================

This is an example of showing how to extract topics from
20 newsgroup dataset with Latent Dirichlet Allocation.
It uses both online and batch update methods to train LDA
model and shows that online method can get meaningful topics
with less iterations.

The code is modified from the
"Topics extraction with Non-Negative Matrix Factorization"
example.

"""

# Authors: Chyi-Kwei Yau

from __future__ import print_function
from time import time
from sklearn.datasets import fetch_20newsgroups
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.decomposition import LatentDirichletAllocation


# training document number
n_samples = 2000
# holdout documents number
n_holdout = 1000

n_features = 1000
n_top_words = 10
n_topics = 10

# Training iterations for online and batch update.
# online method can converge with less trainign iterations.
n_online_iters = 5
n_batch_iters = 10


def print_top_words(model, feature_names, n_top_words):
    for topic_idx, topic in enumerate(model.components_):
        print("Topic #%d:" % topic_idx)
        print(" ".join([feature_names[i] for i in topic.argsort()[:-n_top_words - 1:-1]]))
    print()



print("Loading 20 news groups dataset...")
dataset = fetch_20newsgroups(shuffle=True, random_state=1,
                             remove=('headers', 'footers', 'quotes'))
vectorizer = CountVectorizer(max_df=0.8, max_features=n_features,
                             min_df=3, stop_words='english')

lda_online = LatentDirichletAllocation(n_topics=n_topics, max_iter=n_online_iters,
                                       learning_method='online', learning_offset=100.,
                                       random_state=0)

lda_batch = LatentDirichletAllocation(n_topics=n_topics, max_iter=n_batch_iters,
                                      learning_method='batch', random_state=0)

doc_list = dataset.data[:n_samples]
holdout_doc_list = dataset.data[n_samples: n_samples+n_holdout]

doc_mtx = vectorizer.fit_transform(doc_list)
holdout_doc_mtx = vectorizer.transform(holdout_doc_list)
feature_names = vectorizer.get_feature_names()

print("\nFitting LDA models with online update for %d iterations..." % n_online_iters)
t0 = time()
lda_online.fit(doc_mtx)
print("done in %0.3fs." % (time() - t0))

print("\nCalculating hold out document perplexity...")
perplexity_online = lda_online.perplexity(holdout_doc_mtx)
print("hold out document perplexity: %0.3f" % perplexity_online)

print("\nTopics in online LDA model:")
print_top_words(lda_online, feature_names, n_top_words)

print("\nFitting LDA models with batch update for %d iterations..." % n_batch_iters)
t1 = time()
lda_batch.fit(doc_mtx)
print("done in %0.3fs." % (time() - t1))

print("\nCalculating hold out document perplexity...")
perplexity_batch = lda_batch.perplexity(holdout_doc_mtx)
print("hold out document perplexity: %0.3f" % perplexity_batch)

print("\nTopics in batch LDA model:")
print_top_words(lda_batch, feature_names, n_top_words)
