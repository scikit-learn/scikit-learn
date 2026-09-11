import sys
from time import time

from sklearn.datasets import fetch_20newsgroups
from sklearn.feature_extraction.text import CountVectorizer, HashingVectorizer

categories = [
    "alt.atheism",
    "comp.graphics",
    "comp.sys.ibm.pc.hardware",
    "misc.forsale",
    "rec.autos",
    "sci.space",
    "talk.religion.misc",
]

print("Loading 20 newsgroups training data")
raw_data, _ = fetch_20newsgroups(categories=categories, return_X_y=True)
data_size_mb = sum(len(s.encode("utf-8")) for s in raw_data) / 1e6
print(f"{len(raw_data)} documents - {data_size_mb:.3f}MB")

if len(sys.argv) > 1:
    n_jobs = {"n_jobs": int(sys.argv[1])}
else:
    n_jobs = {}

t0 = time()
vectorizer = CountVectorizer(**n_jobs)
vectorizer.fit_transform(raw_data)
duration = time() - t0
print(f"CountVectorizer done in {duration:.3f} s at {data_size_mb / duration:.1f} MB/s")
print(f"Found {len(vectorizer.get_feature_names_out())} unique terms")

t0 = time()
vectorizer = HashingVectorizer(n_features=2 * 18, **n_jobs)
vectorizer.fit_transform(raw_data)
duration = time() - t0
print(
    f"HashingVectorizer done in {duration:.3f} s at {data_size_mb / duration:.1f} MB/s"
)
