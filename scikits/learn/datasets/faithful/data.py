# encoding: utf-8
"""Old faithful dataset."""

__docformat__ = 'restructuredtext'

COPYRIGHT   = """See SOURCE. Pr Azzalini has explicitely given his consent for
the use of those data in scipy."""
TITLE       = "Old Faithful Geyser Data"
SOURCE      = """AZZALINI A., BOWMAN A. W. (1990).  A look at some data on the
Old Faithful Geyser.  Applied Statistics (Journal of the Royal Statistical
Society series C), vol. 39, pp. 357-365. Data collected by the Yellowstone Park
geologist, R. A. Hutchinson.

References: 
    - Härdle, W. (1991) Smoothing Techniques with Implementation in S.  New
    York: Springer.
    - Azzalini, A. and Bowman, A. W. (1990).  A look at some data on the Old
Faithful geyser. Applied Statistics, 39, 357--365.

Those data are exactly the ones from Azzalini and Bowman's article."""

DESCRSHORT  = """Waiting time between eruptions and the duration of the
eruption for the Old Faithful geyser in Yellowstone National Park, Wyoming,
USA. Waiting times and duration time are in seconds"""

DESCRLONG   = """According to Azzalini and Bowman's article, those data
were recorded continuously from 1th August to 15th August 1985.

Some of the durations times are labelled as L, M or S (Large, Small, Medium).
According to Azzalini and Bowman's paper: 'because the unbroken sequence
required measurements to be taken at night, some duration times are recorded as
L (long), S (short) and M (medium). Other data sets do not contain a con-
tinuous stream of data, making it difficult to deal with time series features'
"""

NOTE = """Eruptions time in minutes, waiting time to next eruption in
minutes"""

import numpy as np

def load():
    """load the actual data and returns them.
    
    :returns:
        waiting: array
             waiting time until next eruption, in minutes
        duration: array
             duration of eruption, in minutes
    """
    from faithful import waiting, duration
    waiting = np.array(waiting, np.float)
    duration = np.array(duration, np.float)
    return waiting, duration
