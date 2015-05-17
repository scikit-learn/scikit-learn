"""
==================================================
Topics extraction with Latent Dirichlet Allocation
==================================================

This is an example of showing how to extract topics from
20 newsgroup dataset with Latent Dirichlet Allocation.
Also, it shows the difference between batch and online udpate.

The code is modified from the
"Topics extraction with Non-Negative Matrix Factorization"
example.

"""

# Authors: Chyi-Kwei Yau

from sklearn.datasets import fetch_20newsgroups
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.utils import gen_batches
from sklearn.decomposition import LatentDirichletAllocation


def lda_example():
    """
    Example for LDA online update
    """

    # chunk_size is how many documents we want to use
    # in each online iteration
    chunk_size = 2000
    n_features = 1000
    n_top_words = 10

    print('Example of LDA with online update')
    print("Loading 20 news groups dataset...")
    dataset = fetch_20newsgroups(shuffle=True, random_state=1,
                                 remove=('headers', 'footers', 'quotes'))
    total_docs = len(dataset.data)

    vectorizer = CountVectorizer(max_df=0.8, max_features=n_features,
                                 min_df=3, stop_words='english')

    lda = LatentDirichletAllocation(learning_decay=0.7, learning_offset=512., n_jobs=-1,
                                    total_samples=total_docs, random_state=0, verbose=0)

    for chunk_no, idx_slice in enumerate(gen_batches(total_docs, chunk_size)):
        doc_list = dataset.data[idx_slice]
        if chunk_no == 0:
            doc_mtx = vectorizer.fit_transform(doc_list)
            feature_names = vectorizer.get_feature_names()
        else:
            doc_mtx = vectorizer.transform(doc_list)

        # fit model
        print("\nFitting LDA models with online update on chunk %d..." %
              chunk_no)
        lda.partial_fit(doc_mtx)
        print("Topics after training chunk %d:" % chunk_no)
        for topic_idx, topic in enumerate(lda.components_):
            print("Topic #%d:" % topic_idx)
            print(" ".join([feature_names[i] for i in topic.argsort()[:-n_top_words - 1:-1]]))


if __name__ == '__main__':
    lda_example()
