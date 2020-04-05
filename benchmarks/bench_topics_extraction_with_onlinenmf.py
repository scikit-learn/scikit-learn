"""
===========================================
Benchmark Non-negative Matrix Factorization
===========================================

This is a benchmark of :class:`sklearn.decomposition.NMF` on a corpus
of documents and extract additive models of the topic structure of the
corpus.  The output is a list of topics, each represented as a list of
terms (weights are not shown).

Non-negative Matrix Factorization is applied with the generalized
Kullback-Leibler divergence equivalent to Probabilistic Latent
Semantic Indexing.

The time complexity is polynomial in NMF.

"""

# Author: Olivier Grisel <olivier.grisel@ensta.org>
#         Lars Buitinck
#         Chyi-Kwei Yau <chyikwei.yau@gmail.com>
#         Chiara Marmo <chiara.marmo@inria.fr>
# License: BSD 3 clause

from time import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec

import zipfile as zp
from bs4 import BeautifulSoup

from sklearn.feature_extraction.text import TfidfVectorizer, CountVectorizer
from sklearn.decomposition.nmf_original import NMFOriginal
from sklearn.decomposition import NMF

n_samples = range(10000, 20000, 2000)
n_features = range(2000, 10000, 2000)
batch_size = range(400, 1000, 200) 
n_components = 10

# Load the The Blog Authorship Corpus dataset
# from http://u.cs.biu.ac.il/~koppel/BlogCorpus.htm
# and vectorize it.

print("Loading dataset...")
t0 = time()
with zp.ZipFile("/home/parietal/cmarmo/bench/blogs.zip") as myzip:
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

fig = plt.figure(constrained_layout=True, figsize=(22, 13))

spec = gridspec.GridSpec(ncols=len(n_features), nrows=len(batch_size),
                         figure=fig)

ylabel = "Convergence time"
xlabel = "n_samples"

ax = []

for bj in range(len(batch_size)):
    miny = 999999
    maxy = 0
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

            # Fit the NMF model with Kullback-Leibler divergence
            print("Fitting the NMF model "
                  "(generalized Kullback-Leibler divergence) "
                  "with tf-idf features, n_samples=%d and n_features=%d..."
                  % (n_samples[i], n_features[j]))
            t0 = time()
            nmf = NMFOriginal(n_components=n_components, random_state=1,
                              beta_loss='kullback-leibler', solver='mu',
                              max_iter=1000, alpha=.1, l1_ratio=.5).fit(tfidf)
            timesKL[i] = time() - t0
            print("done in %0.3fs." % (timesKL[i]))

            # Fit the NMF model KL
            print("Fitting the online NMF model (generalized Kullback-Leibler "
                  "divergence) with "
                  "tf-idf features, n_samples=%d and n_features=%d..."
                  % (n_samples[i], n_features[j]))
            t0 = time()
            minibatch_nmf = NMF(n_components=n_components,
                                batch_size=batch_size[bj],
                                random_state=1, beta_loss='kullback-leibler',
                                solver='mu', max_iter=1000, alpha=.1,
                                l1_ratio=.5).fit(tfidf)
            timesmbKL[i] = time() - t0
            print("done in %0.3fs." % (timesmbKL[i]))

        ax.append(fig.add_subplot(spec[bj,j], xlabel=xlabel, ylabel= ylabel))
        plt.grid(True)

        str1 = "NMF"
        str2 = "Online NMF"
        ax_index = j+bj*(len(n_features)-1)
        ax[ax_index].plot(n_samples, timesKL, marker='o', label = str1)
        ax[ax_index].plot(n_samples, timesmbKL, marker='o', label = str2)

        ax[ax_index].xaxis.set_major_formatter(ticker.EngFormatter())

        strdesc = "n_features " + str(n_features[j])

        miny = min(miny, min(timesKL), min(timesmbKL))
        maxy = max(maxy, max(timesKL), max(timesmbKL))

        ax[ax_index].set_title(strdesc)

    for j in range(len(n_features)):
        ax_index = j+bj*(len(n_features)-1)
        ax[ax_index].set_ylim(miny-10, maxy+10)

    ax[bj*(len(n_features)-1)+1].legend(bbox_to_anchor=(1.05, 1),
                                        loc='upper left', borderaxespad=0.)
    strbatch = "batch size: " + str(batch_size[bj]) + \
               "\nn_components: " + str(n_components)
    ax[bj*(len(n_features)-1)+1].annotate(strbatch, (1.05, 0.5),
                                          xycoords='axes fraction',
                                          va='center')

plt.savefig('bench_topics.png')
#plt.show()
