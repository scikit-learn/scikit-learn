#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Last Change: Mon Jul 09 03:00 PM 2007 J

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

"""Iris dataset."""

__docformat__ = 'restructuredtext'

COPYRIGHT   = """See SOURCE. """
TITLE       = "Heart disease dataset"
SOURCE      = """http://www.liacc.up.pt/ML/old/statlog/datasets.html"""

DESCRSHORT  = """Heart disease of 270 patients with 13 corresponding features."""

DESCRLONG   = DESCRSHORT

NOTE        = """
Number of Instances: 270 (150 for class has heart disease, 12- for no heart disease).

Number of Attributes: 13.

Attribute Information:

    - age
    - sex
    - chest pain type (4 values)
    - resting blood pressure
    - serum cholestoral in mg/dl
    - fasting blood sugar > 120 mg/dl
    - resting electrocardiographic results (values 0,1,2)
    - maximum heart rate achieved
    - exercise induced angina
    - oldpeak = ST depression induced by exercise relative to rest
    - the slope of the peak exercise ST segment
    - number of major vessels (0-3) colored by flourosopy
    - thal: 3 = normal; 6 = fixed defect; 7 = reversable defect 

Missing Attribute Values: None

Attributes types:

    - Real 1,4,5,8,10,12
    - Ordered:11,
    - Binary: 2,6,9
    - Nominal:7,3,13

Class Distribution: 33.3% for each of 3 classes.
"""

def load():
    """load the heart data and returns them.
    
    :returns:
        data: recordarray
            a record array of the data.
    """
    raise NotImplementedError("This dataset is not available yet, sorry")
    import numpy
    descr = [('age', numpy.int), ('sex', numpy.int),
            ('chest pain', numpy.int), ('resting blood pressure', numpy.int),
            ('serum cholesterol', numpy.float), ('high fasting blood sugar', numpy.int),
            ('resting ECG', numpy.int), ('maximum heart rate', numpy.float),
            ('angina', numpy.int), ('oldpeak', numpy.float),
            ('slope', numpy.int), ('major vessels colored', numpy.int),
            ('thal', numpy.int)]
