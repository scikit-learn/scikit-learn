"""
===========================================
Benchmark Non-negative Matrix Factorization
===========================================

This is a benchmark of :class:`sklearn.decomposition.NMF` on a corpus
of documents and extract additive models of the topic structure of the
corpus.  The output is a list of topics, each represented as a list of
terms (weights are not shown).

Non-negative Matrix Factorization is applied with two different objective
functions: the Frobenius norm, and the generalized Kullback-Leibler divergence.
The latter is equivalent to Probabilistic Latent Semantic Indexing.

The default parameters (n_samples / n_features / n_components) should make
the example runnable in a couple of tens of seconds. You can try to
increase the dimensions of the problem, but be aware that the time
complexity is polynomial in NMF.

"""

# Author: Olivier Grisel <olivier.grisel@ensta.org>
#         Lars Buitinck
#         Chyi-Kwei Yau <chyikwei.yau@gmail.com>
# License: BSD 3 clause

from time import time
import numpy as np
import matplotlib.pyplot as plt

from sklearn.feature_extraction.text import TfidfVectorizer, CountVectorizer
from sklearn.decomposition import NMF
from sklearn.datasets import fetch_20newsgroups

n_samples = range(1000, 1000, 1000)
n_features = range(500, 2500, 1000)
batch_size = 1000
n_components = 10
n_top_words = 20


def print_top_words(model, feature_names, n_top_words):
    for topic_idx, topic in enumerate(model.components_):
        message = "Topic #%d: " % topic_idx
        message += " ".join([feature_names[i]
                             for i in topic.argsort()[:-n_top_words - 1:-1]])
        print(message)
    print()


# Load the 20 newsgroups dataset and vectorize it. We use a few heuristics
# to filter out useless terms early on: the posts are stripped of headers,
# footers and quoted replies, and common English words, words occurring in
# only one document or in at least 95% of the documents are removed.

print("Loading dataset...")
t0 = time()
data, _ = fetch_20newsgroups(shuffle=True, random_state=1,
                             remove=('headers', 'footers', 'quotes'),
                             return_X_y=True)
print("done in %0.3fs." % (time() - t0))

ax1 = plt.subplot(221, ylabel = "time")
ax2 = plt.subplot(222, xlabel = "n_samples", ylabel = "time", sharex = ax1)
ax3 = plt.subplot(223, sharex = ax1, sharey = ax1)
ax3 = plt.subplot(224, xlabel = "n_samples", sharex = ax1, sharey = ax1)


for j in range(len(n_features)):
  timesFr = np.zeros(len(n_samples))
  timesmbFr = np.zeros(len(n_samples))
  timesKL = np.zeros(len(n_samples))
  timesmbKL = np.zeros(len(n_samples))

  for i in range(len(n_samples)):
    data_samples = data[:n_samples[i]]
    # Use tf-idf features for NMF.
    print("Extracting tf-idf features for NMF...")
    tfidf_vectorizer = TfidfVectorizer(max_df=0.95, min_df=2,
                                       max_features=n_features[j],
                                       stop_words='english')
    t0 = time()
    tfidf = tfidf_vectorizer.fit_transform(data_samples)
    print("done in %0.3fs." % (time() - t0))

    # Fit the NMF model
    print("Fitting the NMF model (Frobenius norm) with tf-idf features, "
          "n_samples=%d and n_features=%d..."
          % (n_samples[i], n_features[j]))
    t0 = time()
    nmf = NMF(n_components=n_components, random_state=1,
              alpha=.1, l1_ratio=.5).fit(tfidf)
    timesFr[i] = time() - t0
    print("done in %0.3fs." % (timesFr[i]))

    # Fit the NMF model with minibatch
    print("Fitting the online NMF model (Frobenius norm) with tf-idf features, "
          "n_samples=%d and n_features=%d..."
          % (n_samples[i], n_features[j]))
    t0 = time()
    minibatch_nmf = NMF(n_components=n_components, batch_size=batch_size,
                        random_state=1, alpha=.1, l1_ratio=.5,
                        max_iter=3).fit(tfidf)
    timesmbFr[i] = time() - t0
    print("done in %0.3fs." % (timesmbFr[i]))

    # Fit the NMF model
    print("Fitting the NMF model (generalized Kullback-Leibler divergence) with "
          "tf-idf features, n_samples=%d and n_features=%d..."
          % (n_samples[i], n_features[j]))
    t0 = time()
    nmf = NMF(n_components=n_components, random_state=1,
              beta_loss='kullback-leibler', solver='mu', max_iter=1000, alpha=.1,
              l1_ratio=.5).fit(tfidf)
    timesKL[i] = time() - t0
    print("done in %0.3fs." % (timesKL[i]))

    # Fit the NMF model
    print("Fitting the NMF model (generalized Kullback-Leibler divergence) with "
          "tf-idf features, n_samples=%d and n_features=%d..."
          % (n_samples[i], n_features[j]))
    t0 = time()
    minibatch_nmf = NMF(n_components=n_components, batch_size=batch_size,
                       random_state=1, beta_loss='kullback-leibler',
                       solver='mu', max_iter=1000, alpha=.1,
                       l1_ratio=.5).fit(tfidf)
    timesmbKL[i] = time() - t0
    print("done in %0.3fs." % (timesmbKL[i]))

  str1 = "Features " + str(n_features[j]) 
  ax1.plot(n_samples, timesFr)
  ax2.plot(n_samples, timesKL)
  ax3.plot(n_samples, timesmbFr, label = str1 )

ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

plt.subplots_adjust(wspace=0, hspace=0)
plt.show()
