#! /usr/bin/env python
# Last Change: Fri Jun 08 11:00 AM 2007 J

# This script generates a python file from the txt data
import csv

dataname = 'Oldfaithful.txt'
f = open(dataname, 'r')
a = csv.reader(f)
el = [i for i in a]
duration = [i[0] for i in el]
waiting = [i[1] for i in el]

# Convert duration and waiting times in second
duration2 = []
for i in range(len(duration)):
	if duration[i] == 'L':
		duration2.append('L')
	elif duration[i] == 'M':
		duration2.append('M')
	elif duration[i] == 'S':
		duration2.append('S')
	else:
		m, s = duration[i].split(':')
		m = int(m)
		s = int(s)
		assert s >= 0 and s < 60
		duration2.append(m * 60 + s)
waiting2 = [int(i) * 60 for i in waiting]

# Write the data in oldfaitful.py
a = open("oldfaithful.py", "w")

a.write("duration =  [\n")
for i in range(len(duration2) - 1):
	a.write("'%s', " % duration2[i])
a.write("'%s']\n" % duration2[-1])

a.write("waiting =  [\n")
for i in range(len(waiting2) - 1):
	a.write("%s, " % waiting2[i])
a.write("%s]\n" % waiting2[-1])
a.close()
