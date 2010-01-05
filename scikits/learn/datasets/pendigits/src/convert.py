#! /usr/bin/env python
# Last Change: Fri Jun 22 05:00 PM 2007 J

# This script generates a python file from the txt data
import csv

dataname_tra = 'pendigits.tra'
dataname_tes = 'pendigits.tes'
f = open(dataname_tra, 'r')
a = csv.reader(f)
tra = [[int(j) for j in i] for i in a]

f = open(dataname_tes, 'r')
a = csv.reader(f)
tes = [[int(j) for j in i] for i in a]

# Write the data in pendigits.py
a = open("pendigits.py", "w")

a.write("training =  [\n")
for i in range(len(tra)-1):
    a.write(tra[i].__repr__() + ',\n')
a.write(tra[-1].__repr__() + ']\n\n')

a.write("testing =  [\n")
for i in range(len(tes)-1):
    a.write(tes[i].__repr__() + ',\n')
a.write(tes[-1].__repr__() + ']\n')

a.close()
