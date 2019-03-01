from time import time

from scipy import sparse
import pandas as pd

from sklearn.decomposition.nmf import _beta_divergence
from sklearn.feature_extraction.text import CountVectorizer, HashingVectorizer

from nmf import NMF
from nmf_original import NMFOriginal

import matplotlib.pyplot as plt
from dirty_cat.datasets import fetch_traffic_violations

dataset = 'traffic_violations'

try:
    X = sparse.load_npz('X.npz')
except FileNotFoundError:
    if dataset == 'wiki':
        df = pd.read_csv('/home/pcerda/parietal/online_nmf/scikit-learn/' +
                         'enwiki_1000000_first_paragraphs.csv')
        cats = df['0'].astype(str)
        counter = HashingVectorizer(analyzer='word', ngram_range=(1, 1),
                                    n_features=2**12, norm=None,
                                    alternate_sign=False)
    elif dataset == 'traffic_violations':
        data = fetch_traffic_violations()
        df = pd.read_csv(data['path'])
        cats = df['Model'].astype(str).values
        counter = CountVectorizer(analyzer='char', ngram_range=(3, 3))
    X = counter.fit_transform(cats)
    # sparse.save_npz('X.npz', X)

n_test = 10000
n_train = 50000

X_test = X[:n_test, :]
X = X[n_test:n_train + n_test, :]

n_components = 10

print(X.shape)

time_nmf = []
kl_nmf = []
time_nmf2 = []
kl_nmf2 = []

fig, ax = plt.subplots()
# plt.yscale('log')
fontsize = 16
beta_loss = 'kullback-leibler'

max_iter_nmf = [1, 5, 10, 30, 50, 100]
max_iter_minibatch_nmf = [1, 5, 10, 20, 30, 40]

nmf2 = NMF(
    n_components=n_components, beta_loss=beta_loss, batch_size=1000,
    solver='mu', max_iter=1, random_state=10, tol=0)

for i, max_iter in enumerate(zip(max_iter_nmf, max_iter_minibatch_nmf)):
    nmf = NMFOriginal(n_components=n_components, beta_loss=beta_loss,
                      solver='mu', max_iter=max_iter[0], random_state=10,
                      tol=0)
    t0 = time()
    nmf.fit(X)
    W = nmf.transform(X_test)
    tf = time() - t0
    time_nmf.append(tf)
    print('Time NMF: %.1fs.' % tf)
    kldiv = _beta_divergence(X_test, W, nmf.components_,
                             nmf.beta_loss) / X_test.shape[0]
    kl_nmf.append(kldiv)
    print('KL-div NMF: %.2f' % kldiv)
    del W

    t0 = time()
    # nmf2 = NMF(
    #     n_components=n_components, beta_loss=beta_loss, batch_size=1000,
    #     solver='mu', max_iter=max_iter[1], random_state=10, tol=0)
    nmf2.partial_fit(X)
    W = nmf2.transform(X_test)
    tf = time() - t0
    time_nmf2.append(tf)
    print('Time MiniBatchNMF: %.1fs.' % tf)
    kldiv = _beta_divergence(X_test, W, nmf2.components_,
                             nmf2.beta_loss) / X_test.shape[0]
    kl_nmf2.append(kldiv)
    print('KL-div MiniBatchNMF: %.2f' % kldiv)
    del W

    if i > 0:
        plt.plot(time_nmf, kl_nmf, 'r', marker='o')
        plt.plot(time_nmf2, kl_nmf2, 'b', marker='o')
        plt.pause(.01)
        if i == 1:
            plt.legend(labels=['NMF', 'Online NMF'], fontsize=fontsize)


plt.tick_params(axis='both', which='major', labelsize=fontsize-2)
plt.xlabel('Time (seconds)', fontsize=fontsize)
plt.ylabel(beta_loss, fontsize=fontsize)

if dataset == 'traffic_violations':
    title = 'Traffic Violations; Column: Model'
elif dataset == 'wiki':
    title = 'Wikipedia articles (first paragraph)'
ax.set_title(title, fontsize=fontsize+4)

figname = 'benchmark_nmf_%s.pdf' % dataset
print('Saving: ' + figname)
plt.savefig(figname,
            transparent=False, bbox_inches='tight', pad_inches=0)
plt.show()
