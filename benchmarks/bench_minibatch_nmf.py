from time import time

from sklearn.decomposition._nmf import _beta_divergence
from sklearn.utils import gen_batches

import zipfile as zp
from bs4 import BeautifulSoup

from sklearn.feature_extraction.text import TfidfVectorizer

from sklearn.decomposition import NMF, MiniBatchNMF, non_negative_factorization

import matplotlib.pyplot as plt
import matplotlib.lines as mlines


def get_optimal_w(X, H):
    W, _, _ = non_negative_factorization(
        X=X, W=None, H=H,
        n_components=n_components,
        init='custom', update_H=False, solver='mu',
        beta_loss=beta_loss, tol=1e-4, max_iter=200, alpha=0.,
        l1_ratio=0., regularization=None, random_state=None,
        verbose=0, shuffle=False)
    return W


n_components = 10
n_features = 500
beta_loss = 'kullback-leibler'
n_train = 12000
n_test = 7000
batch_sizes = [1000, 2000, 4000]
forget_factors = [1., 0.5]
random_state = 12
color = ['b', 'g', 'c', 'm', 'y', 'k']

# Load the The Blog Authorship Corpus dataset
# from http://u.cs.biu.ac.il/~koppel/BlogCorpus.htm
# and vectorize it.

print("Loading dataset...")
t0 = time()
with zp.ZipFile("/home/cmarmo/software/test/blogs.zip") as myzip:
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

# Use tf-idf features for NMF.
print("Extracting tf-idf features for NMF...")
tfidf_vectorizer = TfidfVectorizer(max_df=0.95, min_df=2,
                                   max_features=n_features,
                                   stop_words='english')
t0 = time()
X = tfidf_vectorizer.fit_transform(data)
print("done in %0.3fs." % (time() - t0))

X_test = X[:n_test, :]
X = X[n_test:n_train + n_test, :]

max_iter_nmf = [1, 5, 10, 30, 50, 100]
n_iter_minibatch_nmf = 50

fig, ax = plt.subplots()
plt.xscale('log')
fontsize = 10

c = 0
labels = []
handles = []

for batch_size in batch_sizes:

    n_batch = (n_train - 1) // batch_size + 1

    for forget_factor in forget_factors:

        minibatch_nmf = MiniBatchNMF(
            n_components=n_components, beta_loss=beta_loss,
            batch_size=batch_size,
            solver='mu', random_state=random_state, max_iter=3,
            forget_factor=forget_factor)

        total_time = 0
        time_nmf = []
        loss_nmf = []

        labels.append(('MiniBatchNMF '
                       f'{batch_size= }'
                       f' {forget_factor= }'))
        handles.append(mlines.Line2D([], [], color=color[c], marker='o'))

        for n_iter in range(n_iter_minibatch_nmf):

            for j, slice in enumerate(
                gen_batches(n=n_train,
                            batch_size=batch_size)
                           ):
                t0 = time()
                minibatch_nmf.partial_fit(X[slice])
                tf = time() - t0
                total_time += tf
                if ((j % 11 == 9) and (n_iter <= 1)) or j == n_batch - 1:
                    time_nmf.append(total_time)
                    W = get_optimal_w(X_test, minibatch_nmf.components_)
                    loss = _beta_divergence(X_test, W,
                                            minibatch_nmf.components_,
                                            minibatch_nmf.beta_loss) / n_test
                    loss_nmf.append(loss)
                    plt.plot(time_nmf, loss_nmf, color=color[c], alpha=0.3,
                             linestyle='-', marker='o',
                             label=labels[-1])
                    plt.pause(.01)

            print('Time MiniBatchNMF: %.1fs.' % total_time)
            print('KL-div MiniBatchNMF: %.2f' % loss)
            del W

        c += 1

total_time = 0
time_nmf = []
loss_nmf = []
for i, max_iter in enumerate(max_iter_nmf):
    nmf = NMF(n_components=n_components, beta_loss=beta_loss,
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
    plt.plot(time_nmf, loss_nmf, 'r', marker='o', label='NMF')
    plt.pause(.01)
    del W

labels.append('NMF')
handles.append(mlines.Line2D([], [], color='r', marker='o'))

plt.legend(handles=handles, labels=labels, fontsize=fontsize-2)
plt.tick_params(axis='both', which='major', labelsize=fontsize-2)
plt.xlabel('Time (seconds)', fontsize=fontsize)
plt.ylabel(beta_loss, fontsize=fontsize)
title = ('Blog Authorship Corpus dataset')
ax.set_title(title, fontsize=fontsize+4)

figname = 'benchmark_nmf_blog_authorship.png'
print('Saving: ' + figname)
plt.savefig(figname, transparent=False)
plt.show()
