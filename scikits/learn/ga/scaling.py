#genetic algorithm scaling routines
#based on galib.
#
from ga_util import *
import scipy.stats as stats
from Numeric import *
from scipy_base.fastumath import *
# if a score is less the 2 standard deviations below, the average, its score
# is arbitrarily set to zero
class sigma_truncation_scaling:
	def __init__(self,scaling = 2):
		self.scaling = scaling	
	def scale(self,pop):
		sc = pop.scores()	
		avg = my_mean(sc)
		if len(sc) > 1: dev = my_std(sc)
		else: dev = 0
		f = sc - avg + self.scaling * dev
		f=choose(less_equal(f,0.),(f,0.))
		for i in range(len(pop)): pop[i].fitness(f[i])
		return pop	

class no_scaling:
	def scale(self,pop): 
		for ind in pop: ind.fitness(ind.score())
		return pop	

class linear_scaling:
	def __init__(self,mult = 1.2):
		self.mult = mult
	def scale(self,pop):
		sc = pop.scores()
		pmin = min(sc)
		if pmin < 0: raise GAError, 'linear scaling does not work with objective scores < 0'
		pmax = max(sc)
		pavg = my_mean(sc)
		if(pavg == pmax): 
			a = 1.
			b = 0.
		elif pmin > (self.mult * pavg - pmax)/(self.mult - 1.):
			delta = pmax - pavg
			a = (self.mult - 1.) * pavg / delta
			b = pavg * (pmax - self.mult * pavg) / delta
		else:
			delta = pavg - pmin
			a = pavg / delta
			b = -pmin * pavg / delta
		f = sc * a + b
		f=choose(less_equal(f,0.),(f,0.))
		for i in range(len(pop)): pop[i].fitness(f[i])
