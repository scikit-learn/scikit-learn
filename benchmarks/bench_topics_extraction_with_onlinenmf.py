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

from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.decomposition.nmf_original import NMFOriginal
from sklearn.decomposition import NMF

n_samples = range(10000, 20000, 2000)
n_features = range(2000, 10000, 2000)
batch_size = 600
n_components = range(10, 70, 20)

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

spec = gridspec.GridSpec(ncols=len(n_features), nrows=len(n_components),
                         figure=fig)

ylabel = "Convergence time"
xlabel = "n_samples"

ax = []

for bj in range(len(n_components)):
    miny = 999999
    maxy = 0
    for j in range(len(n_features)):
        timesKL = np.zeros(len(n_samples))
        timesmbKL = np.zeros(len(n_samples))
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
            nmf = NMFOriginal(n_components=n_components[bj], random_state=1,
                              beta_loss='kullback-leibler', solver='mu',
                              max_iter=1000, alpha=.1, l1_ratio=.5).fit(tfidf)
            timesKL[i] = time() - t0
            print("done in %0.3fs." % (timesKL[i]))
            lossKL[i] = nmf.reconstruction_err_

            # Fit the NMF model KL
            print("Fitting the online NMF model (generalized Kullback-Leibler "
                  "divergence) with "
                  "tf-idf features, n_samples=%d and n_features=%d..."
                  % (n_samples[i], n_features[j]))
            t0 = time()
            minibatch_nmf = NMF(n_components=n_components[bj],
                                batch_size=batch_size,
                                random_state=1, beta_loss='kullback-leibler',
                                solver='mu', max_iter=1000, alpha=.1,
                                l1_ratio=.5).fit(tfidf)
            timesmbKL[i] = time() - t0
            print("done in %0.3fs." % (timesmbKL[i]))
            lossmbKL[i] = minibatch_nmf.reconstruction_err_

        ax.append(fig.add_subplot(spec[bj, j], xlabel=xlabel, ylabel=ylabel))
        plt.grid(True)

        str1 = "time NMF"
        str2 = "time Online NMF"
        str3 = "loss NMF"
        str4 = "loss Online NMF"

        ax_index = j+bj*len(n_features)
        ax[ax_index].plot(n_samples, timesKL, marker='o', label=str1)
        ax[ax_index].plot(n_samples, timesmbKL, marker='o', label=str2)

        ax2 = ax[ax_index].twinx()
        ax2.set_ylabel('loss')

        ax2.plot(n_samples, lossKL, marker='x', ls='dashed', label=str3)
        ax2.plot(n_samples, lossmbKL, marker='x', ls='dashed', label=str4)

        ax[ax_index].xaxis.set_major_formatter(ticker.EngFormatter())

        strdesc = "n_features " + str(n_features[j])

        miny = min(miny, min(timesKL), min(timesmbKL))
        maxy = max(maxy, max(timesKL), max(timesmbKL))

        ax[ax_index].set_title(strdesc)

    for j in range(len(n_features)):
        ax_index = j+bj*len(n_features)
        ax[ax_index].set_ylim(miny-10, maxy+10)

    ax[(bj+1)*len(n_features)-1].legend(bbox_to_anchor=(1.05, 1),
                                        loc='upper left', borderaxespad=0.)
    ax2.legend(bbox_to_anchor=(1.05, 1),
               loc='lower left', borderaxespad=0.)
    strbatch = "batch size: " + str(batch_size) + \
               "\nn_components: " + str(n_components[bj])
    ax[(bj+1)*len(n_features)-1].annotate(strbatch, (1.05, 0.8),
                                          xycoords='axes fraction',
                                          va='center')

plt.savefig('bench_topics.png')
# plt.show()
