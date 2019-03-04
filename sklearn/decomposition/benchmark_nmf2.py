
from time import time

from scipy import sparse
import pandas as pd

from sklearn.decomposition.nmf import _beta_divergence
from sklearn.feature_extraction.text import CountVectorizer, HashingVectorizer
from sklearn.utils import gen_batches

from nmf import NMF
from nmf_original import NMFOriginal
from nmf_original import non_negative_factorization

import matplotlib.pyplot as plt

limit = 1000000
j = 0
articles = []
file = 'enwiki_1M_first_paragraphs.csv'
for i, line in enumerate(open('enwiki_preprocessed_with_articles_markup.txt')):
    if line.startswith('<article'):
        article = ''
        print(line)
        continue
    if line.startswith('</article>'):
        articles.append(article)
        continue
    if article == '':
        article = line
    if len(articles) >= limit:
        break
df = pd.DataFrame(articles)
df.to_csv('%d_first_paragraphs.csv' % len(articles))

# Donload file from:
# https://filesender.renater.fr/?s=download&token=88222d6d-5aee-c59b-4f34-c233b4d184e1
df = pd.read_csv('enwiki_1000000_first_paragraphs.csv')
cats = df['0'].sample(frac=1, random_state=5).astype(str)
counter = HashingVectorizer(analyzer='word', ngram_range=(1, 1),
                            n_features=2**12, norm=None,
                            alternate_sign=False)
X = counter.fit_transform(cats)
n_components = 10
beta_loss = 'kullback-leibler'
n_train = 300000
n_test = 10000
batch_size = 10000
random_state = 12
n_batch = (n_train - 1) // batch_size + 1
X_test = X[:n_test, :]
X = X[n_test:n_train + n_test, :]

max_iter_nmf = [1, 5, 10, 30, 50, 100]
n_iter_minibatch_nmf = 10


def get_optimal_w(X, H):
    W, _, _ = non_negative_factorization(
        X=X, W=None, H=H,
        n_components=n_components,
        init='custom', update_H=False, solver='mu',
        beta_loss=beta_loss, tol=1e-4, max_iter=200, alpha=0.,
        l1_ratio=0., regularization=None, random_state=None,
        verbose=0, shuffle=False)
    return W


minibatch_nmf = NMF(
    n_components=n_components, beta_loss=beta_loss, batch_size=batch_size,
    solver='mu', random_state=random_state, max_iter=3)

fig, ax = plt.subplots()
plt.xscale('log')
fontsize = 16

total_time = 0
time_nmf = []
loss_nmf = []
for n_iter in range(n_iter_minibatch_nmf):

    for j, slice in enumerate(gen_batches(n=n_train,
                                          batch_size=batch_size)):
        t0 = time()
        minibatch_nmf.partial_fit(X[slice])
        tf = time() - t0
        total_time += tf
        if ((j % 11 == 9) and (n_iter == 0)) or j == n_batch - 1:
            time_nmf.append(total_time)
            W = get_optimal_w(X_test, minibatch_nmf.components_)
            loss = _beta_divergence(X_test, W, minibatch_nmf.components_,
                                    minibatch_nmf.beta_loss) / n_test
            loss_nmf.append(loss)
            if j == n_batch - 1:
                plt.plot(time_nmf[-1], loss_nmf[-1],
                         'b', marker='o')
            else:
                plt.plot(time_nmf[-1], loss_nmf[-1],
                         'b', marker='+')
            plt.pause(.01)

    print('Time MiniBatchNMF: %.1fs.' % total_time)
    print('KL-div MiniBatchNMF: %.2f' % loss)
    del W

total_time = 0
time_nmf = []
loss_nmf = []
for i, max_iter in enumerate(max_iter_nmf):
    nmf = NMFOriginal(n_components=n_components, beta_loss=beta_loss,
                      solver='mu', max_iter=max_iter,
                      random_state=random_state, tol=0)
    t0 = time()
    nmf.fit(X)
    tf = time() - t0
    total_time += tf
    time_nmf.append(total_time)
    print('Time NMF: %.1fs.' % total_time)
    W = get_optimal_w(X_test, nmf.components_)
    loss = _beta_divergence(X_test, W, nmf.components_,
                            nmf.beta_loss) / n_test
    loss_nmf.append(loss)
    print('KL-div NMF: %.2f' % loss)
    plt.plot(time_nmf, loss_nmf, 'r', marker='o')
    plt.pause(.01)
    del W

plt.legend(labels=['NMF', 'Mini-batch NMF'], fontsize=fontsize)
plt.tick_params(axis='both', which='major', labelsize=fontsize-2)
plt.xlabel('Time (seconds)', fontsize=fontsize)
plt.ylabel(beta_loss, fontsize=fontsize)
title = 'Wikipedia articles (first paragraph)'
ax.set_title(title, fontsize=fontsize+4)

figname = 'benchmark_nmf_wikipedia_articles.pdf'
print('Saving: ' + figname)
plt.savefig(figname,
            transparent=False, bbox_inches='tight', pad_inches=0)
plt.show()
