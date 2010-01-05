#! /usr/bin/env python
# Last Change: Tue Jul 17 05:00 PM 2007 J

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
a = open("../pendigits.py", "w")

a.writelines(dumpvar(tra, 'training'))
a.writelines(dumpvar(tes, 'testing'))
a.close()
