#base definitions for genetic algorithms
import scipy.rv
import scipy.stats

GAError = 'GA Error'

def nop(x): return x
def flip_coin(p): return (scipy.rv.random() < p)

import whrandom

def flip_coin2(p): return (whrandom.random() < p)
class empty_class: pass

def shallow_clone(item):
	new = empty_class()
	new.__class__ = item.__class__
	new.__dict__.update(item.__dict__)
	return new
#these are exacly correct, but htey prevent problems with -Inf and Inf	
def my_std(s):
#	try:
		a = remove_NaN(s)
		if len(a) > 1: return scipy.stats.stdev(a)
		else: return 0.
#	except: 
#		import pdb
#		pdb.set_trace()
def my_mean(s):
	a = remove_NaN(s)
	if len(a) > 1: return scipy.stats.mean(a)
	else: return 0.
	
def testflip():
	
	import time
	b = time.clock()
	for i in range(10000): a = flip_coin(.5)
	e = time.clock()
	print 'rv_flip',e-b
	b = time.clock()
	for i in range(10000): a = flip_coin2(.5)
	e = time.clock()
	print 'wh_flip',e-b
	from rv import random
	b = time.clock()
	for i in range(10000): 
		a = random() < .5
	e = time.clock()
	print 'rv',e-b
	from whrandom import random
	b = time.clock()
	for i in range(10000): 
		a = random() < .5
	e = time.clock()
	print 'wh',e-b


#This causes a seg fault on windows (sometimes)...
def remove_NaN(z):
	#global INF, NEG_INF
	INF = 1e300**2 		# These lines are the culprits 
	NEG_INF = -1e300**3 # and they seem to interact strangly with pyGrad and pyPlot
	front = 0; back = -1;
	if(hasattr(z,'tolist')): ss = z.tolist()
	else: ss = z
	#if NEG_INF == 0.: NEG_INF = -1e300**3
	#if INF == 0.: INF = 1e300**2
	try: back = ss.index(NEG_INF)
	except: pass
		#import sys
		#print 'error:',sys.exc_type,sys.exc_value, NEG_INF, type(ss)
		#import pdb
		#pdb.set_trace()
	try: front = ss.index(INF)+1
	except ValueError: pass	

	#print front, back, NEG_INF	
	if (front != 0 or back != -1): return ss[front:back]
	else: return ss

