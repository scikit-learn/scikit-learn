#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Last Change: Fri Jun 22 05:00 PM 2007 J

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

"""Pen-based recognition of handwritten digits."""

__docformat__ = 'restructuredtext'

COPYRIGHT   = """See SOURCE. """
TITLE       = """Pen-Based Recognition of Handwritten Digits Original,
unnormalized version"""
SOURCE      = """
E. Alpaydin, Fevzi. Alimoglu
Department of Computer Engineering
Bogazici University, 80815 Istanbul Turkey
alpaydin@boun.edu.tr
September 1998

References: 
	F. Alimoglu (1995) Combining Multiple Classifiers for Pen-Based
	Handwritten Digit Recognition, 
	MSc Thesis, Institute of Graduate Studies in Science and 
	Engineering, Bogazici University.
	http://www.cmpe.boun.edu.tr/~alimoglu/alimoglu.ps.gz

	F. Alimoglu, E. Alpaydin, "Methods of Combining Multiple Classifiers 
	Based on Different Representations for Pen-based Handwriting
	Recognition," Proceedings of the Fifth Turkish Artificial 
	Intelligence and Artificial Neural Networks Symposium (TAINN 96), 
	June 1996, Istanbul, Turkey.
	http://www.cmpe.boun.edu.tr/~alimoglu/tainn96.ps.gz."""

DESCRSHORT  = """
We create a digit database by collecting 250 samples from 44 writers.
The samples written by 30 writers are used for training,
cross-validation and writer dependent testing, and the digits 
written by the other 14 are used for writer independent testing. This
database is available in the UNIPEN format."""

DESCRLONG   = """
We create a digit database by collecting 250 samples from 44 writers.
The samples written by 30 writers are used for training,
cross-validation and writer dependent testing, and the digits 
written by the other 14 are used for writer independent testing. This
database is available in the UNIPEN format.

We use a WACOM PL-100V pressure sensitive tablet with an integrated 
LCD display and a cordless stylus. The input and display areas are
located in the same place. Attached to the serial port of an Intel 
486 based PC, it allows us to collect handwriting samples. The tablet
sends $x$ and $y$ tablet coordinates and pressure level values of the
pen at fixed time intervals (sampling rate) of 100 miliseconds. 

These writers are asked to write 250 digits in random order inside 
boxes of 500 by 500 tablet pixel resolution.  Subject are monitored 
only during the first entry screens. Each screen contains five boxes
with the digits to be written displayed above. Subjects are told to
write only inside these boxes.  If they make a mistake or are unhappy
with their writing, they are instructed to clear the content of a box 
by using an on-screen button. The first ten digits are ignored 
because most writers are not familiar with this type of input devices,
but subjects are not aware of this. 

In our study, we use only ($x, y$) coordinate information. The stylus
pressure level values are ignored. The raw data that we capture 
from the tablet consist of integer values between 0 and 500 
(tablet input box resolution). 

pendigits-orig contain original, unnormalized data. pendigits is the 
normalized and resampled version where all inputs are of the same
length. Here because of speed or the digit, feature vectors may be of
different lengths, e.g., '1' is shorter than '8'.
"""

NOTE        = """
Number of Instances
-------------------

pendigits-orig.tra	Training	7494
pendigits-orig.tes	Testing		3498

The way we used the dataset was to use first half of training for 
actual training, one-fourth for validation and one-fourth
for writer-dependent testing. The test set was used for 
writer-independent testing and is the actual quality measure.

pendigits.tra divides into
    Actual training	set	3748
    Validation		1873
    Writer-dependent test	1873

Number of Attributes
--------------------

Input size depends on writing speed and time and is not fixed +1 class
attribute

For Each Attribute:
-------------------

The data is in the UNIPEN format. See
I. Guyon UNIPEN 1.0 Format Definition, 
ftp://ftp.cis.upenn.edu/pub/UNIPEN-pub/definition/unipen.def
1994

Missing Attribute Values
------------------------

Class Distribution
------------------

          classes
          0    1    2    3    4    5    6    7    8    9
tra     384  390  392  370  391  375  351  375  363  357 Tot 3748 
cv      209  201  201  163  185  163  191  196  178  186 Tot 1873 
wdep    187  188  187  186  204  182  178  207  178  176 Tot 1873 
windep  363  364  364  336  364  335  336  364  336  336 Tot 3498 
"""

def load():
    """load the actual data and returns them.
    
    :returns:
        data: dictionary.
            data['training'] and data['testing'] both are dictionaries which
            contains the following keys: x, y and class. x and y are the
            coordinates of the eight resampled points, and class is an integer
            between 0 and 9, indicating the number label.

    example
    -------

    Let's say you want to plot the first sample of the training set with
    matplotlib. You would do something like plot(data['training']['x'][0],
    data['training']['y'][0], '-')
    """
    import numpy
    from pendigits import training, testing
    assert len(training) == 7494
    assert len(testing) == 3498

    def raw_to_num(dt):
        coordinates = numpy.empty((len(dt), 16), numpy.int)
        digclass = numpy.empty(len(dt), dtype = numpy.int)
        for i in range(len(coordinates)):
            coordinates[i] = dt[i][:-1]
            digclass[i] = dt[i][-1]
        xcor = coordinates[:, ::2]
        ycor = coordinates[:, 1::2]
        return xcor, ycor, digclass

    xcor, ycor, digclass = raw_to_num(training)
    training = {'x' : xcor, 'y' : ycor, 'class' : digclass}
    xcor, ycor, digclass = raw_to_num(testing)
    testing = {'x' : xcor, 'y' : ycor, 'class' : digclass}
    return {'testing' : testing, 'training' : training}
