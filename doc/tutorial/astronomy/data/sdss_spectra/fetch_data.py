import os
import urllib2
import numpy as np

DATA_URL = ('http://www.astro.washington.edu/users/'
            'vanderplas/pydata/spec4000_corrected.npz')
LOCAL_FILE = 'spec4000_corrected.npz'

# data directory is password protected so the public can't access it    
password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()
password_mgr.add_password(None, DATA_URL, 'pydata', 'astroML')
handler = urllib2.HTTPBasicAuthHandler(password_mgr)
opener = urllib2.build_opener(handler)

# download training data
if not os.path.exists(LOCAL_FILE):
    print "downloading data from", DATA_URL
    fhandle = opener.open(DATA_URL)
    open(LOCAL_FILE, 'w').write(fhandle.read())
