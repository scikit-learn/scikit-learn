"""Script to download the movie review dataset"""

from pathlib import Path
from hashlib import sha256
import tarfile
from urllib.request import urlopen


URL = "http://www.cs.cornell.edu/people/pabo/movie-review-data/review_polarity.tar.gz"

ARCHIVE_SHA256 = "fc0dccc2671af5db3c5d8f81f77a1ebfec953ecdd422334062df61ede36b2179"
ARCHIVE_NAME = Path(URL.rsplit("/", 1)[1])
DATA_FOLDER = Path("txt_sentoken")


if not DATA_FOLDER.exists():

    if not ARCHIVE_NAME.exists():
        print("Downloading dataset from %s (3 MB)" % URL)
        opener = urlopen(URL)
        with open(ARCHIVE_NAME, "wb") as archive:
            archive.write(opener.read())

    try:
        print("Checking the integrity of the archive")
        assert sha256(ARCHIVE_NAME.read_bytes()).hexdigest() == ARCHIVE_SHA256

        print("Decompressing %s" % ARCHIVE_NAME)
        with tarfile.open(ARCHIVE_NAME, "r:gz") as archive:
            archive.extractall(path=".")

    finally:
        ARCHIVE_NAME.unlink()
