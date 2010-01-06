#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Last Change: Sun Jul 22 03:00 PM 2007 J

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

In order to train and test our classifiers, we need to represent 
digits as constant length feature vectors. A commonly used technique
leading to good results is resampling the ( x_t, y_t) points. 
Temporal resampling (points regularly spaced in time) or spatial
resampling (points regularly spaced in arc length) can be used here. 
Raw point data are already regularly spaced in time but the distance
between them is variable. Previous research showed that spatial
resampling to obtain a constant number of regularly spaced points 
on the trajectory yields much better performance, because it provides 
a better alignment between points. Our resampling algorithm uses 
simple linear interpolation between pairs of points. The resampled
digits are represented as a sequence of T points ( x_t, y_t )_{t=1}^T,
regularly spaced in arc length, as opposed to the input sequence, 
which is regularly spaced in time.

So, the input vector size is 2*T, two times the number of points
resampled. We considered spatial resampling to T=8,12,16 points in our
experiments and found that T=8 gave the best trade-off between 
accuracy and complexity.
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

16 (8 (x, y) coordinates)

For Each Attribute:
-------------------

All input attributes are integers in the range 0..100.

Missing Attribute Values
------------------------

None 
"""
