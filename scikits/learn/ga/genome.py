"""
**WARNING: SETTING MUTATION THROUGH THE INITIALIZE METHOD IS BROKEN RIGHT NOW.**

The class genome is used as the base class for genome classes that you
use in your program.  It should not be used directly.  In particular, the **clone()** 
method should be overridden in the derived class.  Also companion
classes for mutation, initialization, crossover, and evaluation must be defined 
for derived genomes. 	

The module also contains the class **list_genome**.  This genome is very similar to 
the standard genome used in genetic algorithms.  It consist of a series of genes
concatenated together to form a list.  The mutation, initialization and crossover
operators are already defined for this class.  All you need to do is assign 
it an evaluator class and you are in business.
"""

from ga_util import *
import scipy.stats as rv
import Numeric, copy 
import tree

class default_evaluator: 
	""" This default evaluator class just reminds you to define your own. """
	def evaluate(self,genome): 
		return genome.performance() #if a performance() method is available, use it!
		#raise GAError, 'objective must be specified'

class genome:
	"""
	The class genome is used as the base class for genome classes that you
	use in your program.  It should not be used directly.  In particular, the **clone()**
	and **copy()** classes should be overridden in a derived class.  Also companion
	classes for mutation, initialization, crossover, and evaluation must be defined 
	for derived genomes. 	

	The **score()** method returns the raw score of the genome.  The **fitness()** function 
	returns the scaled score for the genome.  The scaled score must be set by the population
	class object.
	
	evaluated -- boolean flag (0 or 1).  Indicates whether the genome has been evaluated.
		       It is useful to prevent mutliple evaluations of a genome that was evaluated
			 in a previous generation and was not changed in the current generation.
	evals -- integer.  Keeps track of how many times the genome was evaluated *this generation*.
		   The population resets this value to 0 every generation and then takes inventory
		   of how many genomes were evaluated during after the evaluation phase using this
		   parameter.  It is kind of a screwed up way to keep track of the number of evaluations,
		   but it works.
	_score -- number.  The *raw score* for the genome.  It should be accessed through 
		    the **score()** function.
	_fitness --	number.  The *scaled score* for the genome.  It should be accessed through 
		      the **fitness()** function.   
	"""
	evaluator = default_evaluator()
	def __init__(self):
		self.evaluated = 0
		self.evals=0
	def initialize(self,settings = None):
		"""initialize the genome.  This calls the **initializer.evaluate()** 
 		   function to initialize the genome.
		   This method is usually called by the population initializer.
		
		   settings -- This is a dictionary of genetic algorithm parameters
		  		   passed in from the population initializer.  The only
				   setting that is used here is 'p_mutate'.  The value from
				   the 'p_mutate' key is passed to each gene's **set_mutation()**
				   function.  If the 'p_mutate' key is not present, then
				   set_mutation is not called.
		   NOTE: SETTING THE MUTATION RATE IS BROKEN RIGHT NOW		   
		"""
		self.evaluated = 0
		self.evals = 0
		self.initializer.evaluate(self)
	def clone(self): 
		""" **The clone method must be overridden in a base class.** 
		   It should return a new copy of the genome.  Take care that
		   the new genome does not share data with the old genome.  This
		   could result in hard to find bugs.  On the other hand, this
		   method is called frequently.  If it is slow, it will definitely
		   affect the speed of your code, so try to copy things in an
		   efficient way.
		"""
		raise GAError, 'must define clone() in your genome class'
#	def copy(self,other): 
#		""" **The clone method must be overridden in a base class.** """	
#		raise GAError, 'must define copy() in your genome class'
	def touch(self): 
		"""Resets the evaluated flag to 0 so that the genome will be re-evaluated
		   next time **score()** is called.
		"""
		self.evaluated = 0
	def mutate(self):
		"""Calls the **mutator.evaluate()** function to mutate this genome.  Most
		"""
		if self.mutator.evaluate(self):
			self.evaluated = 0
			return 1
		return 0	
	def evaluate(self, force = 0):
		"""If *evaluated = 0* then, the **evaluator.evaluate()** function is
		   called to evaluate the genome.  The evaluated flag is set to 1, 
		   and evals is incremented.
		    		    		
		   Arguments: 

		   force -- boolean flag (0 or 1).  This forces the evaluation of the
		            genome even if *evaluated = 1*.
		"""
		if (not self.evaluated) or force:
			self._score = self.evaluator.evaluate(self)
			self.evaluated = 1
			self.evals = self.evals + 1
		return self._score	
	def score(self,*val):
		""" Returns the current *raw score* of for the genome.  If the genome is not
		    evaluated, it evaluates the genome before returning the score.  It
		    can also be used to set the score for the genome.
		    
		    val -- number (optional).  If val is specified, the raw score for the
		           gene is set to val and evaluated is set to 1.
		"""
		if len(val): 
			self._score = val[0]
			self.evaluated = 1
		else:	self.evaluate()
		return self._score
	def fitness(self,*val):
		""" Returns or set the current *scaled score* for the genome.  
		    
		    val -- number (optional).  If val is specified, the scaled score for the
		           gene is set to val.
		"""
		if len(val): self._fitness = val[0]
		return self._fitness
	def validate(self): 
		"""overload this function to check that the genome
		   has a valid structure.  For example you could
		   reject it if it had more than a certain number of
		   node_types or leaves.  You can also reject things
		   that have more than a certain depth. etc.
		"""
		return 1
		
class list_genome_default_initializer:
	""" The evaluate() function for this class simply calls the **initialize()**
	    function for each gene in the **list_genome**.
	"""
	def evaluate(self,genome): 
		for gene in genome: gene.initialize()
	def __call__(self,genome): return self.evaluate(genome)		
		
class list_genome_default_mutator:		
	""" The evaluate() function for this class simply calls the **mutate()**
	    function for each gene in the **list_genome**. It returns 1 if
	    any of the genes were mutated
	"""
	def evaluate(self,genome): 
		mutated = 0
		for gene in genome: mutated = gene.mutate() or mutated 
		return mutated
	def __call__(self,genome): return self.evaluate(genome)		
		
class list_genome_singlepoint_crossover:
	def evaluate(self,parents):
		#assume mom and dad are the same length
		mom = parents[0]; dad = parents[1]
		if(len(mom) > 1): 
		    crosspoint = rv.randint(1,len(mom)-1).rvs()[0]
		else: 
		    crosspoint = rv.randint(0,len(mom)).rvs()[0]
		brother = (mom[:crosspoint] + dad[crosspoint:]).clone()
		sister = (dad[:crosspoint] + mom[crosspoint:]).clone()
		return brother, sister
	def __call__(self,parents): return self.evaluate(parents)		

import ga_list
class list_genome(genome,ga_list.ga_list):
	""" This genome is very similar to the standard genome used in genetic algorithms.  
	    It consist of a series of genes concatenated together to form a list.  The mutation, 
	    initialization and crossover operators are already defined for this class.  
	    All you need to do is assign it an evaluator class and you are in business.
	    
	    The list of genes is managed by the ga_list base class.  Most of the standard
	    list operations can be used slice, concatenate , and get items from the genome.
	    Care must be taken to when using these operators to call **touch()** after using
	    these operators in a way that alters the genome.  This assures that the _score of
	    the genome coincides with the current gene values of the genome.
	     
	    You can also customize this genome by defining your own initializer, crossover
	    and mutation classes.
	   
	    crossover -- crossover object.  Defaults to a single-point crossover object.
	    mutator -- mutator object.  Default simply calls the genes' mutator.
	    initializer -- initializer object.  Default simply calls the genes' initializer.
	"""    

	default_mutator = list_genome_default_mutator
	default_initializer = list_genome_default_initializer
	singlepoint_crossover = list_genome_singlepoint_crossover

	crossover = singlepoint_crossover()
	mutator = default_mutator()
	initializer = default_initializer()
	def __init__(self,list = None):
		"""Arguments:
			 list -- list of genes (optional).  If this is not specified, you
		               can build the genome using **append()** to add genes.
		"""
		genome.__init__(self)
		ga_list.ga_list.__init__(self,list)
	def initialize(self,settings = None):
		genome.initialize(self,settings)
		if settings and settings.has_key('p_mutate'):
			for g in self: g.set_mutation(settings['p_mutate'])
	def clone(self): 
		"""This returns a new genome object.  The new genome is a shallow copy 
		   of all the object's attributes.  The gene's in the ga_list are all cloned.
		   If your gene has an attribute that is a dictionary or some other complex
		   object, you may need to override this function so that it explicitly copies
		   your complex object.  Otherwise, the clone and the original object will
		   end up sharing data (a bad thing).
		"""
		return ga_list.ga_list.data_clone(self)
	def array(self): 
		"""Most of the time, the genes in this genome specify numeric parameters.
		   This method returns the values of the genes in an array (NumPy) 
		"""
		return Numeric.array(self.get_values())
	def set_values(self,x):
		""" Set the values of the genes
		"""
		for i in range(len(self)):
			self[i].set_value(x[i])
	def get_values(self):
		""" Return the actual vlues of the genes as a list
		"""
		return map(lambda x: x.value(),self)
	"""	
	def pick_numbers(self):
		start = []; lower = []; upper =[];
		for gene in self:
			s = gene._value
			l,u = gene.bounds
			start = start + [s]
			lower = lower + [l]
			upper = upper + [u]
		return start, lower, upper
	"""
def dict_choice(dict):
	tot = 0
	for key in dict.keys(): tot = tot + len(dict[key])
	index = rv.choice(xrange(0,tot))
	for key in dict.keys(): 
		if index >= len(dict[key]):
			index = index - len(dict[key])
		else:
			return key,dict[key][index]	
	#shouldn't get here		
	return None,None

def in_list	(list,val):
	try: 
		list.index(val)
		return 1
	except ValueError: return 0

SymbolError = 'SymbolError'
NoneError = 'NoneError'
class tree_crossover:
	cross_rejects = ['ST']
	def __init__(self):
		self.cross_point = {}
	def bad_cross_point(self,sym):	
		try: 
			self.cross_rejects.index(sym)
			return 1
		except ValueError: return 0		
	def evaluate(self,parents):
		"""Takes a tuple of two parent tree_genomes.  It wraps the
		   function crosser() which does the real crossover work.
		   This function just validates the children produced
		   by crosser(). If they are not valid, the process is 
		   tried over again up to 10 times until it succeeds.
		   If it fails to create valid children from the parents 
		   after 10 tries, a ValueError is raised.  Otherwise
		   the two crossed clones are returned.
		   Note:  If crosser() raises a SymbolError, this is
		   propagated up as a ValueError.
		"""
		for i in range(10):
			try:
				sis,bro = self.crosser(parents)
				if sis.validate() and bro.validate():
					return sis,bro
			except SymbolError: break
			except NoneError: 
				print "hmmm. None for a parent value, try again"
				print "      This often happens when 'ST' isn't included in rejects list"
				
		raise ValueError
	def crosser(self,parents):	   
		"""Takes a tuple of two parent tree_genomes.  Clone the parents.
		   From one parent clone, select a random node making sure 
		   that the nodes derive_type is not in the list of symbols 
		   to reject (cross_rejects).  From the other parent clone,
		   choose another random node THAT HAS THE SAME SYMBOL TYPE 
		   and swap the two nodes.  If the same symbol is not found,
		   a new symbol is chosen from clone1 and the process is tried
		   again up to 10 times.  A SymbolError is raised if the symbol
		   is never found.  Otherwise the two crossed clones are returned.
		"""
		sib1 = parents[0].clone(); sib2 = parents[1].clone()
		sis = sib1.root; bro = sib2.root
		tries = 50 #try 50 times to find a compatible cross symbol
		tried_sym = []
		for i in range(tries):
			sym,node_a = dict_choice(sis.symbol_table)
			if not self.bad_cross_point(sym) and bro.symbol_table.has_key(sym): 
				break
			elif i == (tries - 1): 
				msg = "chosen symbol not found in dad (%s tries)" % `tries`
				raise SymbolError, msg
			else: tried_sym.append(sym)	
		node_b = rv.choice(bro.symbol_table[sym])
		idx = 0
		try:
			for child in node_a.get_parent().children():
				if node_a is child: break
				else: idx = idx + 1		
			node_a.get_parent().children()[idx] = node_b
			idx = 0
			for child in node_b.get_parent().children():
				if node_b is child: break
				else: idx = idx + 1		
		except AttributeError:
			print 'symbol:',sym 		
			raise NoneError
		node_b.get_parent().children()[idx] = node_a
		#now get nodes pointing at the correct parents
		temp = node_a.get_parent()
		node_a.set_parent(node_b.get_parent())
		node_b.set_parent(temp)
		sib1.evaluated = 0; sib2.evaluated = 0
		if self.cross_point.has_key(sym): 
			self.cross_point[sym] =  self.cross_point[sym] + 1
		else: self.cross_point[sym] = 1	
		return sib1,sib2
	def __call__(self,genome): return self.evaluate(genome)	

import language
class tree_genome_default_initializer:
	def evaluate(self,genome): genome.generate()	
	def __call__(self,genome): return self.evaluate(genome)		

class tree_genome_default_mutator:		
	def evaluate(self,genome): return genome.root.mutate()
	def __call__(self,genome): return self.evaluate(genome)	

class tree_genome(genome):
	deleted = 0
	default_mutator = tree_genome_default_mutator
	default_intializer = tree_genome_default_initializer
	default_crossover = tree_crossover
	max_depth = 10000 # just a big number
	mutator = default_mutator()
	initializer = default_intializer()
	crossover = default_crossover()
	def __init__(self,language):
		genome.__init__(self)
		self.root = None
		self.language = language
#	def touch(self): calls genome touch because of inheritance order
	def initialize(self,settings = None):
		genome.initialize(self,settings)
		if settings and settings.has_key('p_mutate'):
			g.root.set_mutation(settings['p_mutate'])
	def defaultize(self):
		""" set the nodes to their default values"""
		if self.root is None:
			genome.initialize(self)
		self.root.defaultize()		
	def generate(self):
		"""Since we have to clean up after circular references, this wraps
		   the language generator and does the clean up.
		"""   
		self.language.max_depth = self.max_depth			
		while 1:
			try: 
				if self.root: self.root.delete_circulars()				
				self.root = self.language.generate()
				self.root.initialize()			
			except language.DepthError: pass
			if self.root and self.validate(): break		
	def clone(self): 
		new = shallow_clone(self)
		if(self.root): new.root = self.root.clone()
		return new		
	def __del__(self):
		"""Because of circular references, I have to specifically
		   delete the root node.
		"""
		#print 'del in'
		if hasattr(self,'root'):
			#print 'del root'
			if self.root: 
				#print 'del circ'
				self.root.delete_circulars()
			del self.root
