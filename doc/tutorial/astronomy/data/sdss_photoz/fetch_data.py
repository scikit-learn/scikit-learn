"""
This file fetches photometric observations associated with SDSS galaxy
spectra which have spectroscopically confirmed redshifts.  This directly
queries the SDSS database for the information, and thus can take a few
minutes to run.
"""

import os
import urllib
import numpy as np

URL = 'http://cas.sdss.org/public/en/tools/search/x_sql.asp'
archive_file = 'sdss_galaxy_colors.npy'

N_objects = 50000

DTYPE = [('u', float),
         ('g', float),
         ('r', float),
         ('i', float),
         ('z', float),
         ('specClass', int),
         ('redshift', float),
         ('redshift_err', float)]

def sql_query(sql_str, url=URL, format='csv'):
    """Execute SQL query"""
    # remove comments from string
    sql_str = ' \n'.join(map(lambda x: x.split('--')[0], sql_str.split('\n')))
    params = urllib.urlencode(dict(cmd=sql_str, format=format))
    return urllib.urlopen(url + '?%s' % params)


query_text = ('\n'.join(
        ("SELECT TOP %i" % N_objects,
         "   p.u, p.g, p.r, p.i, p.z, s.specClass, s.z, s.zerr",
         "FROM PhotoObj AS p",
         "   JOIN SpecObj AS s ON s.bestobjid = p.objid",
         "WHERE ",
         "   p.u BETWEEN 0 AND 19.6",
         "   AND p.g BETWEEN 0 AND 20",
         "   AND s.specClass > 1 -- not UNKNOWN or STAR",
         "   AND s.specClass <> 5 -- not SKY",
         "   AND s.specClass <> 6 -- not STAR_LATE")))


if not os.path.exists(archive_file):
    print "querying for %i objects" % N_objects
    print query_text
    output = sql_query(query_text)
    print "finished.  Processing & saving data"

    data = np.loadtxt(output, delimiter=',', skiprows=1, dtype=DTYPE)
    np.save(archive_file, data)
else:
    print "data already on disk"
