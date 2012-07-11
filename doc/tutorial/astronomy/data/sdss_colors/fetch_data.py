import os
import urllib2
import numpy as np

DTYPE_TRAIN = [('u-g', np.float32),
               ('g-r', np.float32),
               ('r-i', np.float32),
               ('i-z', np.float32),
               ('redshift', np.float32)]

DTYPE_TEST = [('u-g', np.float32),
               ('g-r', np.float32),
               ('r-i', np.float32),
               ('i-z', np.float32),
               ('label', np.int32)]

SDSS_COLORS_URL = "http://www.astro.washington.edu/users/vanderplas/pydata/"
TRAIN_FILE = 'sdssdr6_colors_class_train.dat'
TEST_FILE = 'sdssdr6_colors_class.200000.dat'

# data directory is password protected so the public can't access it    
password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()
password_mgr.add_password(None, SDSS_COLORS_URL, 'pydata', 'astroML')
handler = urllib2.HTTPBasicAuthHandler(password_mgr)
opener = urllib2.build_opener(handler)

# download training data
destination = TRAIN_FILE.rstrip('.dat') + '.npy'
if not os.path.exists(destination):
    url = SDSS_COLORS_URL + TRAIN_FILE
    fhandle = opener.open(url)
    np.save(destination, np.loadtxt(opener.open(url), dtype=DTYPE_TRAIN))

# download test data
destination = TEST_FILE.rstrip('.dat') + '.npy'
if not os.path.exists(destination):
    url = SDSS_COLORS_URL + TEST_FILE
    print "downloading data from", url
    fhandle = opener.open(url)
    np.save(destination, np.loadtxt(opener.open(url), dtype=DTYPE_TEST))

