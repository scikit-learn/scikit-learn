#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Last Change: Sun Jul 01 08:00 PM 2007 J

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
TITLE       = "Iris Plants Database"
SOURCE      = """Creator: R.A. Fisher 
Donor: Michael Marshall (MARSHALL%PLU@io.arc.nasa.gov)
Date: July, 1988

This is a copy of UCI ML iris datasets, except that the data are in mm instead
of cm, so that exact values as int can be given.

References: 
   - Fisher,R.A. "The use of multiple measurements in taxonomic problems"
     Annual Eugenics, 7, Part II, 179-188 (1936); also in "Contributions to
     Mathematical Statistics" (John Wiley, NY, 1950).
   - Duda,R.O., & Hart,P.E. (1973) Pattern Classification and Scene Analysis.
     (Q327.D83) John Wiley & Sons.  ISBN 0-471-22361-1.  See page 218.
   - Dasarathy, B.V. (1980) "Nosing Around the Neighborhood: A New System
     Structure and Classification Rule for Recognition in Partially Exposed
     Environments".  IEEE Transactions on Pattern Analysis and Machine
     Intelligence, Vol. PAMI-2, No. 1, 67-71.
   - Gates, G.W. (1972) "The Reduced Nearest Neighbor Rule".  IEEE Transactions
     on Information Theory, May 1972, 431-433.
   - See also: 1988 MLC Proceedings, 54-64.  Cheeseman et al's AUTOCLASS II
     conceptual clustering system finds 3 classes in the data."""

DESCRSHORT  = """The famous Iris database, first used by Sir R.A Fisher"""

DESCRLONG   = """This is perhaps the best known database to be found in the
pattern recognition literature.  Fisher's paper is a classic in the field and
is referenced frequently to this day.  (See Duda & Hart, for example.)  The
data set contains 3 classes of 50 instances each, where each class refers to a
type of iris plant.  One class is linearly separable from the other 2; the
latter are NOT linearly separable from each other.  """

NOTE        = """
Number of Instances: 150 (50 in each of three classes)

Number of Attributes: 4 numeric, predictive attributes and the class

Attribute Information:
   - sepal length in mm
   - sepal width in mm
   - petal length in mm
   - petal width in mm
   - class: 
        - Iris-Setosa
        - Iris-Versicolour
        - Iris-Virginica

Missing Attribute Values: None

Class Distribution: 33.3% for each of 3 classes.
"""

def load():
    """load the iris data and returns them.
    
    :returns:
        data: recordarray
            a record array of the data.
    """
    import numpy
    from iris import SL, SW, PL, PW, CLI
    PW = numpy.array(PW).astype(numpy.float)
    PL = numpy.array(PL).astype(numpy.float)
    SW = numpy.array(SW).astype(numpy.float)
    SL = numpy.array(SL).astype(numpy.float)
    data    = {}
    for i in CLI.items():
        name = i[0][5:]
        data[name] = numpy.empty(len(i[1]), [('petal width', numpy.int),\
                        ('petal length', numpy.int),
                        ('sepal width', numpy.int),
                        ('sepal length', numpy.int)])
        data[name]['petal width'] = numpy.round(PW[i[1]] * 10)
        data[name]['petal length'] = numpy.round(PL[i[1]] * 10)
        data[name]['sepal width'] = numpy.round(SW[i[1]] * 10)
        data[name]['sepal length'] = numpy.round(SL[i[1]] * 10)
    
    return data
