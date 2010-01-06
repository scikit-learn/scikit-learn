#! /usr/bin/env python
# Last Change: Sun Jul 22 01:00 PM 2007 J

# This script generates a python file from the txt data
import csv
from scikits.learn.datasets.misc import dumpvar

dataname_tra = 'pendigits.tra'
dataname_tes = 'pendigits.tes'
f = open(dataname_tra, 'r')
a = csv.reader(f)
tra = [[int(j) for j in i] for i in a]

f = open(dataname_tes, 'r')
a = csv.reader(f)
tes = [[int(j) for j in i] for i in a]

# Write the data in pendigits.py
ftra = open("../pendigits_tra.py", "w")
ftes = open("../pendigits_tes.py", "w")

ftra.writelines(dumpvar(tra, 'training'))
ftra.close()
ftes.writelines(dumpvar(tes, 'testing'))
ftes.close()
