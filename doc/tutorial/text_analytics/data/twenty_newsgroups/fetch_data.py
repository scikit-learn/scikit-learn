"""Script to download the 20 newsgroups text classification set"""

from pathlib import Path
from hashlib import sha256
import tarfile
from contextlib import closing
from urllib.request import urlopen

URL = "http://people.csail.mit.edu/jrennie/20Newsgroups/20news-bydate.tar.gz"

ARCHIVE_SHA256 = "8f1b2514ca22a5ade8fbb9cfa5727df95fa587f4c87b786e15c759fa66d95610"
ARCHIVE_NAME = Path(URL.rsplit("/", 1)[1])
TRAIN_FOLDER = Path("20news-bydate-train")
TEST_FOLDER = Path("20news-bydate-test")


if not TRAIN_FOLDER.exists() or not TEST_FOLDER.exists():
    if not ARCHIVE_NAME.exists():
        print("Downloading dataset from %s (14 MB)" % URL)
        opener = urlopen(URL)
        with open(ARCHIVE_NAME, "wb") as archive:
            archive.write(opener.read())

    try:
        print("Checking the integrity of the archive")
        assert sha256(ARCHIVE_NAME.read_bytes()).hexdigest() == ARCHIVE_SHA256

        print("Decompressing %s" % ARCHIVE_NAME)
        with closing(tarfile.open(ARCHIVE_NAME, "r:gz")) as archive:
            archive.extractall(path=".")

    finally:
        ARCHIVE_NAME.unlink()
