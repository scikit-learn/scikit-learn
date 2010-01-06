#! /usr/bin/env python
# Last Change: Mon Aug 20 08:00 PM 2007 J

import os
import numpy as N

from scipy.io.arff import loadarff

filename = '../arff.bak/data/pharynx.arff'

cdir = os.path.abspath(os.path.curdir)
wekacmd = "CLASSPATH=" + cdir + os.sep + "weka.jar:." + " java testarff "
a, b = os.popen2(wekacmd + filename)
nattr = int(b.readline())
ninst = int(b.readline())
print "%d nattr and %d instances" % (nattr, ninst)

def parse_nominal(line):
    l = line.split(',', 1)
    name = l[0]
    value = l[1]
    end = value.find('}')
    l = value[1:end].split(',')
    return name, l

def parse_numeric(line):
    l = line.split(',')
    return l[0], float(l[1]), float(l[2]), float(l[3]), float(l[4])

def is_nominal(line):
    if line.find("{") == -1:
        return False
    else:
        return True

haha = []
for line in b.readlines():
    if is_nominal(line):
        haha.append(parse_nominal(line))
    else:
        haha.append(parse_numeric(line))

data, meta = loadarff(filename)

for i in meta:
    print i
    f = data[i]
    print N.min(f)
for i in haha:
    print i
