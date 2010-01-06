#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Last Change: Sun Jul 22 02:00 PM 2007 J

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

"""Pen-based recognition of handwritten digits: testing data"""

__docformat__ = 'restructuredtext'

from scikits.learn.datasets.pendigits import _pendigit as PD

COPYRIGHT   = PD.COPYRIGHT
TITLE       = """Pen-Based Recognition of Handwritten Digits Original: testing data."""
SOURCE      = PD.SOURCE

DESCRSHORT  = PD.DESCRSHORT + """

Those are the testing, unnormalized data."""

DESCRLONG   = PD.DESCRLONG 
NOTE        = PD.NOTE

def load():
    """load the PENDIGITS testing data.
    
    :returns:
        data: dictionary.

    """
    import numpy
    from pendigits_tes import testing
    assert len(testing) == 3498

    def raw_to_num(dt):
        cor = numpy.empty((len(dt), 16), numpy.int)
        digclass = numpy.empty(len(dt), dtype = numpy.int)
        for i in range(len(cor)):
            cor[i] = dt[i][:-1]
            digclass[i] = dt[i][-1]
        #xcor = coordinates[:, ::2]
        #ycor = coordinates[:, 1::2]
        return cor, digclass

    data = {}

    # Coordinates
    cor, digclass = raw_to_num(testing)
    dt = [('x%d' % i, numpy.int) for i in range(8)]
    dt.extend([('y%d' % i, numpy.int) for i in range(8)])
    data['data'] = numpy.empty(len(cor), dt)
    for i in range(8):
        data['data']['x%d' % i] = cor[:, 2*i]
        data['data']['y%d' % i] = cor[:, 2*i + 1]

    # label
    data['label'] = digclass

    # class names
    cnames = ['zero', 'one', 'two', 'three', 'four', 'five', 'six', 'seven',\
              'eight', 'nine']
    data['class'] = numpy.empty(10, 'S%d' % numpy.max([len(i) for i in cnames]))
    for i in range(len(cnames)):
        data['class'][i] = cnames[i]

    return data
