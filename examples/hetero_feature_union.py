"""
============================================
Feature Union with Heterogeneos Data Sources
============================================

Demonstrate combined use of `sklearn.feature_extraction.FeatureUnion` and
`sklearn.pipeline.Pipeline` on a data source with different types of features.
In this example we use the 20-newsgroups dataset and compute features for the
subject line and body in separate pipelines then combine them with a
FeatureUnion and train a learner on the combined set of features.

The choice of features is not particularly helpful, but serves to illustrate
the technique. One could imagine that instead of considering subject & body,
we are considering image and caption or some other heterogeneous datasource.
"""

# Author: Matt Terry <matt.terry@gmail.com>
#
# License: BSD 3 clause

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.datasets import fetch_20newsgroups
from sklearn.datasets.twenty_newsgroups import strip_newsgroup_footer
from sklearn.datasets.twenty_newsgroups import strip_newsgroup_quoting
from sklearn.decomposition import TruncatedSVD
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics import classification_report
from sklearn.pipeline import FeatureUnion
from sklearn.pipeline import Pipeline
from sklearn.svm import SVC


class DictSelector(BaseEstimator, TransformerMixin):
    """Transform a mappable object into a the value at specified key.

    If your data is organized as a dict of sequences, this is helpful for
    selecting data for a specified key.  It is also useful for selecting
    a single column from a pandas DataFrame.

    DictSelects expect a dict of sequences.
    i.e.
    >> data = {'a': [1, 5, 2, 5, 2, 8],
               'b': [9, 4, 1, 4, 1, 3]}
    >> ds = DictSelector(key='a')
    >> data['a'] == ds.transform(data)

    DictSelector expects a dict of sequences.  If you have a sequence of dicts,
    you should use `sklearn.feature_extraction.DictVectorizer`.

    Parameters
    ----------
    key : hashable, required
        The key corresponding to the desired value in a mappable.
    """
    def __init__(self, key):
        self.key = key

    def fit(self, x, y=None):
        return self

    def transform(self, data_dict):
        return data_dict[self.key]


class SubjectBodyExtractor(BaseEstimator, TransformerMixin):
    """Extract the subject & body from a usenet post in a single pass.

    Takes a sequence of strings and produces a dict of sequences.  Keys are
    `subject` and `body`.
    """
    def fit(self, x, y=None):
        return self

    def transform(self, posts):
        subject = []
        body = []
        for text in posts:
            headers, _, bod = text.partition('\n\n')
            bod = strip_newsgroup_footer(bod)
            bod = strip_newsgroup_quoting(bod)
            body.append(bod)

            prefix = 'Subject:'
            sub = ''
            for line in headers.split('\n'):
                if line.startswith(prefix):
                    sub = line[len(prefix):]
                    break
            subject.append(sub)

        return {'subject': subject, 'body': body}


pipeline = Pipeline([
    # Extract the subject & body
    ('subjectbody', SubjectBodyExtractor()),

    # Use FeatureUnion to combine the features from subject and body
    ('union', FeatureUnion([

        # Pipeline for pulling features from the post's subject line
        ('subject', Pipeline([
            ('selector', DictSelector(key='subject')),
            ('tfidf', TfidfVectorizer(min_df=50)),
        ])),

        # Pipeline for pulling features from post's body
        ('body', Pipeline([
            ('selector', DictSelector(key='body')),
            ('tfidf', TfidfVectorizer()),
            ('best', TruncatedSVD(n_components=50)),
        ])),
    ])),

    # Use a SVC classifier on the combined features
    ('svc', SVC(kernel='linear')),
])

train = fetch_20newsgroups(random_state=1, subset='train')
test = fetch_20newsgroups(random_state=1, subset='test')

pipeline.fit(train.data, train.target)
y = pipeline.predict(test.data)
print classification_report(y, test.target)
