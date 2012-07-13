"""
This file fetches photometric observations associated with SDSS galaxy
spectra which have spectroscopically confirmed redshifts.  This directly
queries the SDSS database for the information, and thus can take a few
minutes to run.
"""

import os
import urllib, urllib2
import numpy as np

# Here's how the data can be downloaded directly from the SDSS server.
# This route is limited to N = 50000, so we've done this separately
def fetch_data_sql(N = 50000):
    URL = 'http://cas.sdss.org/public/en/tools/search/x_sql.asp'
    archive_file = 'sdss_galaxy_colors.npy'

    dtype = [('mags', '5float32'),
             ('specClass', 'int8'),
             ('z', 'float32'),
             ('zerr', 'float32')]

    def sql_query(sql_str, url=URL, format='csv'):
        """Execute SQL query"""
        # remove comments from string
        sql_str = ' \n'.join(map(lambda x: x.split('--')[0],
                                 sql_str.split('\n')))
        params = urllib.urlencode(dict(cmd=sql_str, format=format))
        return urllib.urlopen(url + '?%s' % params)

    query_text = ('\n'.join(
            ("SELECT TOP %i" % N,
             "   modelMag_u, modelMag_g, modelMag_r, modelMag_i, modelMag_z, specClass, z, zErr",
             "FROM SpecPhoto",
             "WHERE ",
             "   modelMag_u BETWEEN 0 AND 19.6",
             "   AND modelMag_g BETWEEN 0 AND 20",
             "   AND zerr BETWEEN 0 and 0.03",
             "   AND specClass > 1 -- not UNKNOWN or STAR",
             "   AND specClass <> 5 -- not SKY",
             "   AND specClass <> 6 -- not STAR_LATE")))


    if not os.path.exists(archive_file):
        print "querying for %i objects" % N
        print query_text
        output = sql_query(query_text)
        print "finished.  Processing & saving data"
        try:
            data = np.loadtxt(output, delimiter=',', skiprows=1, dtype=DTYPE)
        except:
            raise ValueError(output.read())
        np.save(archive_file, data)
    else:
        print "data already on disk"


DATA_URL = ('http://www.astro.washington.edu/users/'
            'vanderplas/pydata/sdss_photoz.npy')
LOCAL_FILE = 'sdss_photoz.npy'

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
