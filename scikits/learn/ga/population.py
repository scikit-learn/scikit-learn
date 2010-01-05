#genetic algorithm population
#based on galib. 
import time
import ga_list 
import re, copy
import scaling
import selection
import Numeric
import scipy.stats as stats
from ga_util import *
import pdb

def ftn_minimize(x,y): 
	"""Minimization comparator for fitness (scaled score)."""
	return cmp(x.fitness(),y.fitness())
def ftn_maximize(x,y): 
	"""Maximization comparator for fitness (scaled score)."""
	return cmp(y.fitness(),x.fitness())
def sc_minimize(x,y): 
	"""Minimization comparator for raw score."""
#	return cmp(x.score(),y.score())
	#removed one function call
	return cmp(x.evaluate(),y.evaluate())
def sc_maximize(x,y): 
	"""Maximization comparator for raw score."""
#	return cmp(y.score(),x.score())
	return cmp(y.evaluate(),x.evaluate())
class default_pop_evaluator:
	"""The **evaluate()** method simply calls the **evaluate()**
	   method for all the genomes in the population
	"""
	def evaluate(self,pop,force = 0):
		try:
			evals = 0
			if not pop.evaluated or force:
				for ind in pop: ind.evaluate(force)
		except:
			#this makes where a pop evaluator can simply evaluate a list
			#of genomes - might be useful to simplify remote evaluation
			for ind in pop: ind.evaluate(force)
			
class default_pop_initializer:
	"""The **evaluate()** method simply calls the **evaluate()**
	   method for all the genomes in the population
	"""
	def evaluate(self,pop,settings):
		for i in pop: i.initialize(settings)			
		
class population(ga_list.ga_list):
	"""A population of genomes.  The population is constructed by cloning a 
	   single genome multiple times.
	   
	   The population can be treated, for the most part, like a python list 
	   because it is derived from ga_list.ga_list.  This allows the addition 
	   and subtraction of individuals from the population using list 
	   operators such as slice, append, etc.
	   
         Some of the methods, such as **best()** and **worst()**, assume that 
         the population is sorted.
          
	   **Attributes:**
	   
	   stats -- dictionary of dictionaries. The statistics for the population.
	            This keeps up with such things as the number of evaluations,
	            the best current score, the best ever score, etc.
	   scaler -- An instance of a scaler class.  Used to scale the raw scores of
	             the population to generate the genome fitness values.  If you don't
	             want scaling then use the class **scaling.no_scaling**.  The
	             default scaler is **sigma_truncation_scaling** because it handles
	             negative fitness scores effectively.
	   evaluator -- An instance of an evaluator class.  Used to evaluate the individuals
	                of the population.
	   selector -- An instance of a selector class.  Used to select the individuals
	               of the population that go into the breeding pool.  The default
	               selector is **srs_selector**.
	               
	   **Flags (0 or 1):**
	   
	   evaluated -- Population has been avaluated
	   scaled -- Population has been scaled
	   sorted -- Population is sorted from best to worst
	   select_ready -- The Selector has been updated and is ready for selections
 	   stated -- The population statistics are up to date.            
	"""
	default_scaler = scaling.sigma_truncation_scaling
	default_selector = selection.srs_selector
	default_evaluator =  default_pop_evaluator
	default_initializer =  default_pop_initializer
	
	scaler = default_scaler()
	evaluator = default_evaluator()
	initializer = default_initializer()
#	selector = default_selector()

	def __init__(self,genome,size=1):
		"""Arguments:
		
		   genome -- a genome object.
		   size -- number.  The population size.  The genome will be 
		           replicated size times to fill the population.
		"""
		self.model_genome = genome
		ga_list.ga_list.__init__(self)
		self.ftn_comparator = ftn_maximize
		self.sc_comparator = sc_maximize		
		self._size(size)
		self.selector = population.default_selector() #why'd I do this?
		self.stats={}
	def initialize(self,settings = None):
		"""This method **must** be called before a genetic algorithm 
		   begins evolving the population.  It takes care of initializing
		   the individual genomes, evaluating them, and scaling the population.
		   It also clears and intializes the statistics for the population.
		   
		   Arguments:
		
		   settings -- dictionary of genetic algorithm parameters.  These
		               are passed on to the genomes for initialization.
		"""	
		self.stats = {'current':{},'initial':{},'overall':{}}
		self.stats['ind_evals'] = 0

		print "beigninning genome generation" 
		b = time.clock()	
		self.initializer.evaluate(self,settings)
		e = time.clock()	
		print "finished generation: ", e-b	
		self.touch(); 
		b = time.clock()	
		self.evaluate()
		e = time.clock()	
		print "evaluation time: ", e-b	
		self.scale()
		self.update_stats()
		self.stats['initial']['avg'] = self.stats['current']['avg']
		self.stats['initial']['max'] = self.stats['current']['max']
		self.stats['initial']['min'] = self.stats['current']['min']
		self.stats['initial']['dev'] = self.stats['current']['dev']
	def clone(self): 
		"""Returns a population that has a shallow copy the all the 
		   attributes and clone of all the genomes in the original 
		   object.  It also makes a deep copy of the stats dictionary.
		"""	
		new = ga_list.ga_list.data_clone(self)
		new.stats = {}
		new.stats.update(self.stats)
		return new
	def touch(self):
		"""Reset all the flags for the population."""
		self.evaluated = 0; self.scaled = 0; self.sorted = 0; self.select_ready = 0
		self.stated = 0
	def _size(self, l):
		"""Resize the population."""
		del self[l:len(self)]
		for i in range(len(self),l): self.append(self.model_genome.clone())
		return len(self)		
	def evaluate(self, force = 0):
		"""Call the **evaluator.evaluate()** method to evaluate
		   the population.  The population is also sorted so that 
		   it maintains the correct order.  The population is only 
		   updated if *evaluated=0*.
		   
		   Arguments:
		   
		   force -- forces evaluation even if evaluated = 1
		"""
		b = time.clock()
		self.evaluator.evaluate(self,force)
		e1 = time.clock()
		self.sort()
		e2 = time.clock()
		#this is a cluge to get eval count to work correctly
		preval = self.stats['ind_evals']
		for ind in self:  
			self.stats['ind_evals'] = self.stats['ind_evals'] + ind.evals
			ind.evals = 0		
		print 'evals: ', self.stats['ind_evals'] - preval
		e3 = time.clock()
		self.touch()
		self.evaluated = 1
		e4 = time.clock()
		#print 'eval:',e1-b, 'sort:',e2-e1, 'stats:',e3-e2, 'touch:',e4-e3
	def mutate(self):
		mutations = 0
		for ind in self: mutations  =  mutations + ind.mutate()
		return mutations
	
	def sort(self,type = 'raw', force = 0):
		"""Sort the population so they are ordered from best
		   to worst.  This ordering is specified by the comparator
		   operator used to sort the population.  The comparator
		   is specified usign the **min_or_max()** function. 
		   
		   Arguments:
		   
		   type -- 'raw' or 'scaled'.  Determines wether the
		           sorting is done based on raw scores or on
		           fitness (scaled) scores.
		   force -- forces the sort even if sorted = 1
		"""
#		if not self.sorted or force: 
		if(type == 'scaled'): self.data.sort(self.ftn_comparator)
		elif(type == 'raw'): self.data.sort(self.sc_comparator)	
		else: raise GAError, 'sort type must be "scaled" or "raw"'			
		self.sorted = 1
	def select(self, cnt = 1):
		"""Calls the selector and returns *cnt* individuals.
		
		   Arguments:
		   
		   cnt -- The number of individuals to return.
		"""
		if not self.select_ready:
			self.selector.update(self)
			self.select_ready = 1
		return self.selector.select(self,cnt)
	def scale(self, force = 0):
		"""Calls the **scaler.scale()** method and updates
		   the fitness of each individual.

 		   Arguments:
		   
		   force -- forces the scaling even if scaled = 1
		"""		   
		if not (self.scaled or force):
			self.scaler.scale(self)			
		self.scaled = 1
	def fitnesses(self): 
		"""Returns the fitness (scaled score) of all the
		   individuals in a population as a Numeric array.
		"""		   
		return Numeric.array(map(lambda x: x.fitness(),self))	
	def scores(self):	
		"""Returns the scores (raw) of all the
		   individuals in a population as a Numeric array.
		"""		   
		return Numeric.array(map(lambda x: x.score(),self))	
	def best(self, ith_best = 1): 
		"""Returns the best individual in the population.
		   *It assumes the population has been sorted.*

 		   Arguments:
		   
		   ith_best -- Useful if you want the second(2), third(3), etc.
		               best individual in the population.
		"""		   
		return self[ith_best - 1]
	def worst(self,ith_worst = 1): 
		"""Returns the worst individual in the population.
		   *It assumes the population has been sorted.*

 		   Arguments:
		   
		   ith_worst -- Useful if you want the second(2), third(3), etc.
		               worst individual in the population.
		"""		   
		return self[-ith_worst]		
	def min_or_max(self,*which_one):
		"""Returns or set 'min' or 'max' indicating whether the
		   population is to be minimized or maximized.  
		   *Minimization may require some special handling 
		   in the scaling and selector routines.*

 		   Arguments:
		   
		   which_one -- 'min' or 'max'(optional). Tells the population
		                the problem is a minimization or maximizization
		                problem.
		"""		   
		if len(which_one): 
			if (re.match('min.*',which_one[0],re.I)):
				self.ftn_comparator = ftn_minimize
				self.sc_comparator = sc_minimize				
			elif (re.match('max.*',which_one[0],re.I)):
				self.ftn_comparator = ftn_maximize
				self.sc_comparator = sc_maximize								
			else:	raise GaError, "min_or_max expects 'min' or 'max'"
		if self.ftn_comparator == ftn_minimize: return 'min'
		elif self.ftn_comparator == ftn_maximize: return 'max'
	def update_stats(self):
		"""Update the statistics for the population."""
		s = self.scores()
		self.stats['current']['max'] = max(s)
		self.stats['current']['avg'] = my_mean(s)
		self.stats['current']['min'] = min(s)
		if len(s) > 1: self.stats['current']['dev'] = my_std(s)
		else: self.stats['current']['dev'] = 0	
		try: self.stats['overall']['max'] = max(self.stats['overall']['max'],
								    self.stats['current']['max'])
		except KeyError: self.stats['overall']['max'] = self.stats['current']['max']
		try: self.stats['overall']['min'] = min(self.stats['overall']['min'],
								    self.stats['current']['min'])
		except KeyError: self.stats['overall']['min'] = self.stats['current']['min']
