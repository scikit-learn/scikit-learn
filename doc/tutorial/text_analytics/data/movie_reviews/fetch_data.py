"""Script to download the movie review dataset"""

import os
import tarfile
try:
    from urllib import urlopen
except ImportError:
    from urllib.request import urlopen


URL = ("http://www.cs.cornell.edu/people/pabo/"
       "movie-review-data/review_polarity.tar.gz")

ARCHIVE_NAME = URL.rsplit('/', 1)[1]
DATA_FOLDER = "txt_sentoken"


if not os.path.exists(DATA_FOLDER):

    if not os.path.exists(ARCHIVE_NAME):
        print("Downloading dataset from %s (3 MB)" % URL)
        opener = urlopen(URL)
        with open(ARCHIVE_NAME, 'wb') as archive:
            archive.write(opener.read())

    print("Decompressing %s" % ARCHIVE_NAME)
    try:
        with tarfile.open(ARCHIVE_NAME, "r:gz") as archive:
            archive.extractall(path = '.')
    except AttributeError: # In Python 2.6, tarfile did not yet implement the context manager protocol
        archive = tarfile.open(ARCHIVE_NAME, "r:gz")
        archive.extractall(path = '.')
        archive.close()
    os.remove(ARCHIVE_NAME)
