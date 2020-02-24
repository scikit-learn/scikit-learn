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

import zipfile as zp
from bs4 import BeautifulSoup

from sklearn.feature_extraction.text import TfidfVectorizer, CountVectorizer
from sklearn.decomposition import NMF

n_samples = range(500, 2500, 1000)
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


# Load the The Blog Authorship Corpus dataset
# from http://u.cs.biu.ac.il/~koppel/BlogCorpus.htm
# and vectorize it.

print("Loading dataset...")
t0 = time()
with zp.ZipFile("/home/cmarmo/software/tests/minibatchNMF/blogs.zip") as myzip:
    info = myzip.infolist()
    data = []
    for zipfile in info:
        if not (zipfile.is_dir()):
            filename = zipfile.filename
            myzip.extract(filename)
            with open(filename, encoding='LATIN-1') as fp:
                soup = BeautifulSoup(fp, "lxml")
                text = ""
                for post in soup.descendants:
                    if post.name == "post":
                        text += post.contents[0].strip("\n").strip("\t")
            data.append(text)
print("done in %0.3fs." % (time() - t0))

fig = plt.figure()

ax1 = fig.add_subplot(221, ylabel = "time - gen. KL divergence",
                  title = "standard NMF")
ax1.tick_params(labelbottom=False)
ax2 = fig.add_subplot(222, sharey = ax1,
                  title = "online NMF")
ax2.tick_params(labelbottom=False, labelleft=False)
#ax3 = fig.add_subplot(223, ylabel = "time - Frobenius norm",
#                  xlabel = "n_samples", sharex = ax1)
#ax4 = fig.add_subplot(224, xlabel = "n_samples", sharex = ax2, sharey = ax3)
#ax4.tick_params(labelleft=False)

for j in range(len(n_features)):
  timesFr = np.zeros(len(n_samples))
  timesmbFr = np.zeros(len(n_samples))
  timesKL = np.zeros(len(n_samples))
  timesmbKL = np.zeros(len(n_samples))
  lossFr = np.zeros(len(n_samples))
  lossmbFr = np.zeros(len(n_samples))
  lossKL = np.zeros(len(n_samples))
  lossmbKL = np.zeros(len(n_samples))

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

    # Fit the NMF model Frobenius norm
    #print("Fitting the NMF model (Frobenius norm) with tf-idf features, "
    #      "n_samples=%d and n_features=%d..."
    #      % (n_samples[i], n_features[j]))
    #t0 = time()
    #nmf = NMF(n_components=n_components, random_state=1,
    #          alpha=.1, l1_ratio=.5).fit(tfidf)
    #timesFr[i] = time() - t0
    #print("done in %0.3fs." % (timesFr[i]))

    #print("\nTopics in NMF model:")
    #tfidf_feature_names = tfidf_vectorizer.get_feature_names()
    #print_top_words(nmf, tfidf_feature_names, n_top_words)

    # Fit the NMF model with minibatch Frobenius norm
    #print("Fitting the online NMF model (Frobenius norm) with tf-idf features, "
    #      "n_samples=%d and n_features=%d..."
    #      % (n_samples[i], n_features[j]))
    #t0 = time()
    #minibatch_nmf = NMF(n_components=n_components, batch_size=batch_size,
    #                    random_state=1, alpha=.1, l1_ratio=.5,
    #                    max_iter=3).fit(tfidf)
    #timesmbFr[i] = time() - t0
    #print("done in %0.3fs." % (timesmbFr[i]))

    #print("\nTopics in NMF model:")
    #tfidf_feature_names = tfidf_vectorizer.get_feature_names()
    #print_top_words(nmf, tfidf_feature_names, n_top_words)

    # Fit the NMF model KL
    print("Fitting the NMF model (generalized Kullback-Leibler divergence) with "
          "tf-idf features, n_samples=%d and n_features=%d..."
          % (n_samples[i], n_features[j]))
    t0 = time()
    nmf = NMF(n_components=n_components, random_state=1,
              beta_loss='kullback-leibler', solver='mu', max_iter=1000,
              alpha=.1, l1_ratio=.5).fit(tfidf)
    timesKL[i] = time() - t0
    print("done in %0.3fs." % (timesKL[i]))

    print("\nTopics in NMF model:")
    tfidf_feature_names = tfidf_vectorizer.get_feature_names()
    print_top_words(nmf, tfidf_feature_names, n_top_words)

    # Fit the NMF model KL
    print("Fitting the online NMF model (generalized Kullback-Leibler "
          "divergence) with "
          "tf-idf features, n_samples=%d and n_features=%d..."
          % (n_samples[i], n_features[j]))
    t0 = time()
    minibatch_nmf = NMF(n_components=n_components, batch_size=batch_size,
                       random_state=1, beta_loss='kullback-leibler',
                       solver='mu', max_iter=1000, alpha=.1,
                       l1_ratio=.5).fit(tfidf)
    timesmbKL[i] = time() - t0
    print("done in %0.3fs." % (timesmbKL[i]))

    print("\nTopics in NMF model:")
    tfidf_feature_names = tfidf_vectorizer.get_feature_names()
    print_top_words(nmf, tfidf_feature_names, n_top_words)

  str1 = "n_Ftrs " + str(n_features[j]) 
  ax1.plot(n_samples, timesKL)
  ax2.plot(n_samples, timesmbKL, label = str1)
#  ax3.plot(n_samples, timesFr)
#  ax4.plot(n_samples, timesmbFr)

ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

plt.subplots_adjust(wspace=0, hspace=0, right=0.7)
plt.savefig('bench_topics.png')
plt.show()
