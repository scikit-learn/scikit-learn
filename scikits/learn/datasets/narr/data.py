"""
XXX: currently disabled
"""

"""
=============
NCEP NARR data sets
=============

This module provides access to the National Center for Environment Prediction 
(NCEP) North American Regional Reanalysis (NARR) data sets through the 
Data Access Protocol (DAP). The documentation describing the data set
variables is included in the doc directory. 

See the `pydap home page`_ for information on the DAP client. 

.. _pydap home page: http://pydap.org/
"""

__docformat__ = 'restructuredtext'

TITLE       = "NCEP North American Regional Reanalysis (NARR)"

COPYRIGHT   = """The information on government servers are in the public domain,
unless specifically annotated otherwise, and may be used freely by the public so
long as you do not 1) claim it is your own (e.g. by claiming copyright for NWS 
information -- see below), 2) use it in a manner that implies an endorsement 
or affiliation with NOAA/NWS, or 3) modify it in content and then present it as 
official government material. You also cannot present information of your own 
in a way that makes it appear to be official government information..

The user assumes the entire risk related to its use of this data. NWS is 
providing this data "as is," and NWS disclaims any and all warranties, whether
express or implied, including (without limitation) any implied warranties of 
merchantability or fitness for a particular purpose. In no event will NWS be 
liable to you or to any third party for any direct, indirect, incidental, 
consequential, special or exemplary damages or lost profit resulting from any 
use or misuse of this data.

As required by 17 U.S.C. 403, third parties producing copyrighted works 
consisting predominantly of the material appearing in NWS Web pages must 
provide notice with such work(s) identifying the NWS material incorporated 
and stating that such material is not subject to copyright protection."""

DISCLAIMER = """While every effort has been made to ensure that data are
accurate and reliable within the limits of the current state of
the art, NOAA cannot assume liability for any damages caused by
any errors or omissions in the data, nor as a result of the
failure of the data to function on a particular system.
NOAA makes no warranty, expressed or implied, nor does the fact
of distribution constitute such a warranty.
The data used to produce these analyses and forecasts has
undergone automated quality checks.
"""
 
SOURCE      = """National Centers for Environmental Prediction (NCEP) _`North American Regional Reanalysis (NARR)` _`NOMADS web interface`.

.. __North American Regional Reanalysis (NARR): http://wwwt.emc.ncep.noaa.gov/mmb/rreanl/index.html
.. _NOMADS web interface: http://nomads.ncdc.noaa.gov/data.php#narr_datasets

References: 
    Mesinger, F., et al, 2004: NCEP North American Regional
Reanalysis, 15th Symp. On Global Change and Climate Variations,
Seattle, WA, 11-15 Jan 2004.

    Shafran, P., J. Woollen, W. Ebisuzaki, W. Shi, Y. Fan, R. W.
Grumbine, M. Fennessy, 2004: Observational Data Used for Assimilation in 
the NCEP North American Regional Reanalysis, 20th Intl. Conf. On Interactive
Information Processing Systems for Meteor. Ocean. And Hydrology. Seattle, WA, 11-15
Jan 2004.

    Ebisuzaki, W., J. Alpert, J. Wang, D. Jovic, P. Shafran, 2004: North American
Regional Reanalysis: end user access to large data sets, 20th Intl. Conf. On
Interactive Information Processing Systems for Meteor. Ocean. And Hydrology. Seattle,
WA, 11-15 Jan 2004.

    Mesinger, F., G. DiMego, E. Kalnay, P. Shafran, W. Ebisuzaki, Y. Fan, R. Grumbine,
W. Higgins, Y. Lin, K. Mitchell, D. Parrish, E. Rogers, W. Shi, D. Stokes, J. Woolen,
2003: NCEP Regional Reanalysis, Symp. on Observing and Understanding the Variability
of Water in Weather and Climate, Long Beach, CA, Feb.9-13, 2003.
"""

DESCRSHORT  = """NCEP NARR data sets"""

DESCRLONG   = """The NARR is a reanalysis of
historical observations using a 32 km version of the NCEP 1993
operational ETA model and ETA data assimilation system (EDAS).
The domain of analyses includes North and Central America as
well as small parts of the UK, Eastern Asia and South America
and the oceans in between. The period of the reanalyses is from
October 1978 to the present and analyses were made 8 times
daily. Horizontal boundary conditions were derived from the
NCEP/DOE Reanalysis.
The "merged" dataset provides a high spatial (32 km) and
temporal (3 hour) analyses of North America and adjacent oceans
and land masses from October 1978 to the present. Advantages
over the widely used NCEP/NCAR Reanalysis are its higher
resolution and a much better treatment of the land surface
through a better land-surface model (NOAH), through the
assimilation of more surface data (observed precipitation and
surface winds) and through a better representation of the
terrain (heights, vegetation, soil type).
This data set contains "conventional" atmospheric analyses as
well as model-derived fields which contain estimates of
subsurface, surface, and radiative properties.
"""

import dap.client
import datetime
import numpy as np

NARR_A_URL = "http://nomads.ncdc.noaa.gov:9091/dods/NCEP_NARR_DAILY"

def load(year, month, day, verbose=True):
    """
    Load the National Center for Environment Prediction (NCEP) North American 
    Regional Reanalysis (NARR) data at the given date. The data is available 
    through Data Access Protocol (DAP) from a remote server. The returned 
    object is not the entire data set, but only the metadata describing the 
    data set content. Individual variables are properties of the object, and
    are downloaded only on explicit request (see example below).  
    
    Parameters:
        year : int
        month : int 
        day : int
        
    Note:
    The data is available starting on 1979-01-01 and up to the the month before
    last. Each dataset contains around 190 different variables at 3h intervals. 
    The DAP protocol was formerly known as the Distributed Oceanographic Data 
    Systems (DODS).
    
    
    Some variables included in the set:
        ugrd10m : Wind speed (u) at 10 m altitude [m/s].
        vgrd10m : Wind speed (v) at 10 m altitude [m/s].
        mslet: Mean sea level pressure (eta model) [pa].
        tmpsfc : Surface temperature [K].
        dswrfsfc : Surface downward shortwave radiation flux [W/m^2]. 
        time : Days since 1-1-1  00:00:00.
        lat : Latitude (0 to 89.625) at 0.375 degree resolution [deg].
        lon : Longitude (-220 to -0.625) at 0.375 degree resolution [deg]. 
    
    Example:
        narr = load(2005, 6, 23)
        temp = narr.tmpsfc[:][:][:] # Slices refer to time, lat, lon.
        time = narr.time[:]
        lat = narr.lat[:]
        lon = narr.lon[:]
    
    """

    tail = 'narr-a_221_' 
    dir1 = '%d%02d'%(year, month)
    dir2 = dir1+'%02d'%day
    tail = tail +dir2+'_0000_000'
    url = '/'.join([NARR_A_URL, dir1, dir2, tail])
    if verbose:
        print 'Accessing DODS server: ', url
    return dap.client.open(url)



def convert_time(time):
    """Return datetime objects. 
    
    The NARR dataset provides time in days since 1-1-1 00:00:00:0.0. This 
    function converts the NARR time into a datetime.datetime object. 
    
    Note:
    There seems to be a discrepancy between the days counting of NARR and 
    python datetime. To achieve the fit, two days must be substracted from the
    NARR time.
    """
    time = time-2
    t0 = datetime.datetime(1,1, 1, 0,0,0)
    hour = datetime.timedelta(hours=1)
    day = datetime.timedelta(days=1)
    time = np.atleast_1d(time)
    days = time.astype(int)
    hours =  ((time%1)*24).astype(int)
    return t0+days*day+hours*hour
    
    
# narr = load(1979, 1,1) XXX: currently disabled
