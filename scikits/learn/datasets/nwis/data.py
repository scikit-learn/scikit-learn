#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Last Change: Wed Sep 19 2007

# The code and descriptive text is copyrighted and offered under the terms of
# the BSD License from the authors; see below. However, the actual dataset may
# have a different origin and intellectual property status. See the SOURCE and
# COPYRIGHT variables for this information.

# Copyright (c) 2007 David Huard <david.huard@gmail.com>
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

"""USGS daily discharge at West Branch Delaware River at Walton NY."""

__docformat__ = 'restructuredtext'

COPYRIGHT   = """See src/01423000.dly """
TITLE       = "Daily discharge at West Branch Delaware River at Walton, NY"
SOURCE      = """Creator: National Water Information System, U.S. Geological Survey""" 

DESCRSHORT  = """USGS daily discharge at West Branch Delaware River at Walton, NY."""

DESCRLONG   = """Daily discharge (cubic feet per second) from 2006-9-18 to 2007-9-18 acquired by the 
U.S. Geological Survey and released through the National Water Information System. The data code is included, describing the data quality: 'A' for approved and 'P' for provisional data. Missing discharges, or disharges affected by ice are indicated by the value -999.0. 
"""

NOTE        = """Discharge is computed by measuring the water level and estimating the related discharge through a rating curve. The rating curve is drawn using simultaneously measurements of water level and discharge. The rating curve typically changes over time due to erosion, presence of debris or ice. In this dataset, some discharges are missing due to the presence of ice."""

def load():
    """Return the nwis discharge data.
    
    :Returns:
        d : dict
            contains the following values:
            - 'data' : a record array with the actual data
    
    Example
    -------
    Mask the missing data in discharge and plot the time series.
    
    >>> import datetime, pylab, numpy
    >>> data = load()['data']
    >>> discharge = numpy.ma.masked_values(data['discharge'], -999.0)
    >>> datelist = zip(data['year'], data['month'], data['day'])
    >>> date = [datetime.date(y,m,d) for y,m,d in datelist]
    >>> pylab.plot_date(date, discharge, '-')
    """
    import numpy
    from nwis import date, discharge, code
    import re
    
    # Deal with missing data in discharge
    discharge = numpy.array(discharge)
    missing = (discharge == 'Ice') + (discharge == '')
    discharge[missing] = -999.0
    discharge = discharge.astype(numpy.float)    
    
    # Convert date string to ints
    pattern = r'(\d{4})-(\d{2})-(\d{2})'
    date = numpy.array([re.search(pattern, d).groups(1) for d in date]).astype(numpy.int)
    year, month, day = date.transpose()
        
    data    = {}
    format = [('year', numpy.int), ('month', numpy.int), ('day', numpy.int),
              ('discharge', numpy.float), ('code', str)]
    data['data'] = numpy.empty(len(discharge), format)
    data['data']['year'] = year
    data['data']['month'] = month
    data['data']['day'] = day
    data['data']['discharge'] = discharge   
    return data
