#! /usr/bin/env python
# -*- coding: utf-8 -*-

# The code and descriptive text is copyrighted and offered under the terms of
# the BSD License from the authors; see below. However, the actual dataset may
# have a different origin and intellectual property status. See the SOURCE and
# COPYRIGHT variables for this information.

# Copyright (c) 2007 David Cournapeau <cournape@gmail.com>
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in
#       the documentation and/or other materials provided with the
#       distribution.
#     * Neither the author nor the names of any contributors may be used
#       to endorse or promote products derived from this software without
#       specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Last Change: Fri Jun 08 11:00 AM 2007 J

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
According to Azzalini and Bowman's paper: "because the unbroken sequence
required measurements to be taken at night, some duration times are recorded as
L (long), S (short) and M (medium). Other data sets do not contain a con-
tinuous stream of data, making it difficult to deal with time series features."
"""

NOTE        = """Eruptions time in seconds, waiting time to next eruption (in
seconds)"""

def load():
    """load the actual data and returns them.
    
    :returns:
        data: recordarray
            a record array of the data.
    """
    import numpy
    from oldfaithful import waiting, duration
    assert len(waiting) == len(duration) == 299
    data    = numpy.empty(len(waiting), \
            [('duration', '|S5'), ('waiting', 'int')])
    data['waiting']    = waiting
    data['duration']   = duration
    return data
