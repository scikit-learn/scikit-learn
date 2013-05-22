"""Simple script to fetch a numpy version of the LFW data

Original dataset and credits available at:

  http://vis-www.cs.umass.edu/lfw/

"""
import os
import urllib2

URL = "https://downloads.sourceforge.net/project/scikit-learn/data/lfw_preprocessed.tar.gz"
ARCHIVE_NAME = "lfw_preprocessed.tar.gz"
FOLDER_NAME = "lfw_preprocessed"

if not os.path.exists(FOLDER_NAME):
    if not os.path.exists(ARCHIVE_NAME):
        print "Downloading data, please Wait (58.8MB)..."
        print URL
        opener = urllib2.urlopen(URL)
        open(ARCHIVE_NAME, 'wb').write(opener.read())
        print

    import tarfile
    print "Decompressiong the archive: " + ARCHIVE_NAME
    tarfile.open(ARCHIVE_NAME, "r:gz").extractall()
    os.remove(ARCHIVE_NAME)

