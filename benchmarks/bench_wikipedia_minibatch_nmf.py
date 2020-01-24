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
"""

# Author: Olivier Grisel <olivier.grisel@ensta.org>
#         Lars Buitinck
#         Chyi-Kwei Yau <chyikwei.yau@gmail.com>
# License: BSD 3 clause

from bz2 import BZ2File
import os

from time import time
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt

from scipy import sparse

from urllib.request import urlopen
from joblib import Memory
from sklearn.decomposition import NMF

n_samples = range(1000, 1001, 1)
batch_size = 100
n_components = range(7, 10, 1)

# #############################################################################
# Where to download the data, if not already on disk
redirects_url = "http://downloads.dbpedia.org/3.5.1/en/redirects_en.nt.bz2"
redirects_filename = redirects_url.rsplit("/", 1)[1]

page_links_url = "http://downloads.dbpedia.org/3.5.1/en/page_links_en.nt.bz2"
page_links_filename = page_links_url.rsplit("/", 1)[1]

resources = [
    (redirects_url, redirects_filename),
    (page_links_url, page_links_filename),
]

for url, filename in resources:
    if not os.path.exists(filename):
        print("Downloading data from '%s', please wait..." % url)
        opener = urlopen(url)
        open(filename, 'wb').write(opener.read())
        print()


# #############################################################################
# Loading the redirect files

memory = Memory(location=".")


def index(redirects, index_map, k):
    """Find the index of an article name after redirect resolution"""
    k = redirects.get(k, k)
    return index_map.setdefault(k, len(index_map))


DBPEDIA_RESOURCE_PREFIX_LEN = len("http://dbpedia.org/resource/")
SHORTNAME_SLICE = slice(DBPEDIA_RESOURCE_PREFIX_LEN + 1, -1)


def short_name(nt_uri):
    """Remove the < and > URI markers and the common URI prefix"""
    return nt_uri[SHORTNAME_SLICE]


def get_redirects(redirects_filename):
    """Parse the redirections and build a transitively closed map out of it"""
    redirects = {}
    print("Parsing the NT redirect file")
    for l, line in enumerate(BZ2File(redirects_filename)):
        split = line.split()
        if len(split) != 4:
            print("ignoring malformed line: " + line)
            continue
        redirects[short_name(split[0])] = short_name(split[2])
        if l % 1000000 == 0:
            print("[%s] line: %08d" % (datetime.now().isoformat(), l))

    # compute the transitive closure
    print("Computing the transitive closure of the redirect relation")
    for l, source in enumerate(redirects.keys()):
        transitive_target = None
        target = redirects[source]
        seen = {source}
        while True:
            transitive_target = target
            target = redirects.get(target)
            if target is None or target in seen:
                break
            seen.add(target)
        redirects[source] = transitive_target
        if l % 1000000 == 0:
            print("[%s] line: %08d" % (datetime.now().isoformat(), l))

    return redirects


# disabling joblib as the pickling of large dicts seems much too slow
#@memory.cache
def get_adjacency_matrix(redirects_filename, page_links_filename, limit=None):
    """Extract the adjacency graph as a scipy sparse matrix

    Redirects are resolved first.

    Returns X, the scipy sparse adjacency matrix, redirects as python
    dict from article names to article names and index_map a python dict
    from article names to python int (article indexes).
    """

    print("Computing the redirect map")
    redirects = get_redirects(redirects_filename)

    print("Computing the integer index map")
    index_map = dict()
    links = list()
    for l, line in enumerate(BZ2File(page_links_filename)):
        split = line.split()
        if len(split) != 4:
            print("ignoring malformed line: " + line)
            continue
        i = index(redirects, index_map, short_name(split[0]))
        j = index(redirects, index_map, short_name(split[2]))
        links.append((i, j))
        if l % 1000000 == 0:
            print("[%s] line: %08d" % (datetime.now().isoformat(), l))

        if limit is not None and l >= limit - 1:
            break

    print("Computing the adjacency matrix")
    X = sparse.lil_matrix((len(index_map), len(index_map)), dtype=np.float32)
    for i, j in links:
        X[i, j] = 1.0
    del links
    print("Converting to CSR representation")
    X = X.tocsr()
    print("CSR conversion done")
    return X, redirects, index_map


# stop after 5M links to make it possible to work in RAM
X, redirects, index_map = get_adjacency_matrix(
    redirects_filename, page_links_filename, limit=5000000)
names = {i: name for name, i in index_map.items()}

print(X.shape)

fig = plt.figure()

ax1 = plt.subplot(221)#, ylabel = "time")
ax2 = plt.subplot(222)#, xlabel = "n_samples", ylabel = "time", sharex = ax1)
ax3 = plt.subplot(223)#, sharex = ax1, sharey = ax1)
ax4 = plt.subplot(224)#, xlabel = "n_samples", sharex = ax1, sharey = ax1)


for j in range(len(n_components)):
  timesFr = np.zeros(len(n_samples))
  timesmbFr = np.zeros(len(n_samples))
  timesKL = np.zeros(len(n_samples))
  timesmbKL = np.zeros(len(n_samples))

  for i in range(len(n_samples)):
    X_samples = X[:n_samples[i],:n_samples[i]]

    # Fit the NMF model
    print("Fitting the NMF model (Frobenius norm) on "
          "n_samples=%d and n_components=%d..."
          % (n_samples[i], n_components[j]))
    t0 = time()
    nmf = NMF(n_components=n_components[j], random_state=1,
              alpha=.1, l1_ratio=.5).fit(X_samples)
    timesFr[i] = time() - t0
    print("done in %0.3fs." % (timesFr[i]))

    # Fit the NMF model with minibatch
    print("Fitting the online NMF model (Frobenius norm) on "
          "n_samples=%d and n_components=%d..."
          % (n_samples[i], n_components[j]))
    t0 = time()
    minibatch_nmf = NMF(n_components=n_components[j], batch_size=batch_size,
                        random_state=1, alpha=.1, l1_ratio=.5,
                        max_iter=3).fit(X_samples)
    timesmbFr[i] = time() - t0
    print("done in %0.3fs." % (timesmbFr[i]))

    # Fit the NMF model
    print("Fitting the NMF model (generalized Kullback-Leibler divergence) on "
          "n_samples=%d and n_components=%d..."
          % (n_samples[i], n_components[j]))
    t0 = time()
    nmf = NMF(n_components=n_components[j], random_state=1,
              beta_loss='kullback-leibler', solver='mu', max_iter=1000, alpha=.1,
              l1_ratio=.5).fit(X_samples)
    timesKL[i] = time() - t0
    print("done in %0.3fs." % (timesKL[i]))

    # Fit the NMF model
    print("Fitting the online NMF model (generalized Kullback-Leibler divergence) on "
          "n_samples=%d and n_components=%d..."
          % (n_samples[i], n_components[j]))
    t0 = time()
    minibatch_nmf = NMF(n_components=n_components[j], batch_size=batch_size,
                       random_state=1, beta_loss='kullback-leibler',
                       solver='mu', max_iter=1000, alpha=.1,
                       l1_ratio=.5).fit(X_samples)
    timesmbKL[i] = time() - t0
    print("done in %0.3fs." % (timesmbKL[i]))

  str1 = str(n_components[j]) + " Components"  
  ax1.plot(n_samples, timesFr)
  ax2.plot(n_samples, timesKL)
  ax3.plot(n_samples, timesmbFr, label = str1 )
  ax4.plot(n_samples, timesmbKL)

ax4.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

plt.subplots_adjust(wspace=0, hspace=0)
plt.show()
#fig.savefig('plot.png')
