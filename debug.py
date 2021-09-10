from sklearn.feature_extraction.text import CountVectorizer

corpus = [
    "This is the first document.",
    "This document is the second document.",
    "And this is the third one.",
    "Is this the first document?",
]

vectorizer = CountVectorizer(analyzer="word", min_df=2, max_df=3)
X = vectorizer.fit_transform(corpus)
print(vectorizer.get_feature_names())
# result ['document', 'first']

vectorizer = CountVectorizer(analyzer="word", min_df=2, max_df=3.0)
X = vectorizer.fit_transform(corpus)
print(vectorizer.get_feature_names())
# result ['document', 'first', 'is', 'the', 'this']
