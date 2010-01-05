from ga_util import *
import scipy.stats as stats
rv = stats
#import scipy.io.dumb_shelve
import string, pdb
import os, sys, string
import time, pprint, types,copy
import anydbm, dumbdbm
#import thread, sync
if sys.platform != 'win32': 
	import fcntl
	timer = time.clock	#clock behaves differently work on linux
else:
    timer = time.time
    
dberror = dumbdbm.error

def max_score(pop): return max(map(lambda x: x.score(),pop))

class galg:
	"""A basic genetic algorithm.  The genetic algorithm is responsible
	   for evolving a population of genomes.  While the population and
	   the genomes are in charge of defining most of the genetic operators
	   such as selection, scaling, mutation, and crossover, it is the
	   genetic algorithm class that orchestrates the evolution and calls
	   the operators in the correct order.  Most of the work is done
	   in the **step()** method.
	"""
	valid_settings = ['pop_size','p_replace',
				'p_cross', 'p_mutate','p_deviation',
				'gens','rand_seed','rand_alg','dbase','update_rate']
	output_settings = ['crossover','selector', 'scaler','genome_type']
	default_settings = {'pop_size':150,'p_replace':.8,
				'p_cross': .8, 'p_mutate':'gene',
				'p_deviation': 0.,'gens':35,
				'rand_seed':0,'rand_alg':'CMRG',
				'update_rate': 10000,'dbase':''}
	default_verbose = 1

	def __init__(self,pop):
		self.verbose = self.default_verbose
		self.settings = copy.copy(galg.default_settings)
		self.pop = pop
	def test_settings(self,settings):
		for key in settings.keys(): 
			try: 
				self.output_settings.index(key)
				print 'Warning: The key "%s" in settings is readonly.' % key
			except ValueError:	
				try: self.valid_settings.index(key)
				except ValueError: 
					print 'Warning: The key "%s" in not a valid setting.' % key
					print 'The valid settings are %s' % self.valid_settings
	
	def initialize(self,reseed = 1): 
		b = timer()
		self.test_settings(self.settings)
		self.gen = 0
		sd = self.settings['rand_seed']; alg = self.settings['rand_alg']
		if reseed: rv.initialize(seed = sd, algorithm = alg)
		self.settings['seed_used'] = rv.initial_seed()
		self._print('initializing... seed = %d' % self.settings['seed_used'])
		self.crossover = self.pop.model_genome.crossover # get the crossover op from the first genome
		self.pop.settings = self.settings #should these be shared?
		self.size_pop(self.settings['pop_size'])
			
		self.settings['crossover'] = string.split(str(self.crossover))[0][1:]
		self.settings['selector'] = string.split(str(self.pop.selector))[0][1:]
		self.settings['scaler'] = string.split(str(self.pop.scaler))[0][1:]
		self.settings['genome_type'] = string.split(str(self.pop.model_genome))[0][1:]
#		self._print(self.settings)
		
		self.pop.initialize(self.settings); 
		self.stats = {'selections':0,'crossovers':0,'mutations':0,
				  'replacements':0,'pop_evals':1,'ind_evals':0}
		self.stats.update(self.pop.stats)
		self.step_time = timer() - b
		self.init_dbase()
	def size_pop(self,s):
		self.settings['pop_size'] = s
		self.pop._size(s)

	def step(self,steps=1):
		sz = len(self.pop)
		replace = int(self.settings['p_replace'] * len(self.pop))
		p_crossover = self.settings['p_cross']
		for st in range(steps):
			b = timer()
			for i in range(0,replace,2):
				mom,dad= self.pop.select(2)
				self.stats['selections'] = self.stats['selections'] + 2
				if flip_coin(p_crossover): 
					try: 
						bro,sis = self.crossover((mom,dad))
						self.stats['crossovers'] = self.stats['crossovers'] + 2
						self.pop.append(bro); self.pop.append(sis)				
					except ValueError: 
						#crossover failed - just act as if this iteration never happened
						i = i - 2 
						#print 'crossover failure - ignoring and continuing'
				else: 
					self.pop.append(mom.clone());self.pop.append(dad.clone());
			if replace % 2: #we did one to many - remove the last individual
				del self.pop[-1]
				self.stats['crossovers'] = self.stats['crossovers'] - 1
			e1 = timer();
			self.stats['mutations'] = self.stats['mutations'] + self.pop[sz:].mutate()
#			for ind in self.pop[sz:]:
#				m = ind.mutate()
#				self.stats['mutations'] = self.stats['mutations'] + m
			e2 = timer();
			self.pop.touch()
			self.pop.evaluate()
			e3 = timer();
			del self.pop[sz:] #touch removed from del
			self.pop.scale()
			self.pop.update_stats()
			self.stats['pop_evals'] = self.stats['pop_evals'] + 1
			self.gen = self.gen + 1
			e = timer(); self.step_time = e - b
			#print 'cross:',e1-b,'mutate:',e2-e1,'eval:',e3-e2,'rest:',e-e3
		self.stats.update(self.pop.stats)	
		self.db_entry['best_scores'].append(self.stats['current']['max'])

	def evolve(self):
		b = timer()
		self.initialize()
		self.pre_evolve()
		self.p_dev = self.pop_deviation()
		self.iteration_output()
		while (	self.gen < self.settings['gens'] and
				self.settings['p_deviation'] < self.p_dev  ):
			self.step()
			self.p_dev = self.pop_deviation()
			self.iteration_output()
			if(self.gen % self.settings['update_rate'] == 0):
				self.update_dbase()
		self.update_dbase() #enter status prior to post_evolve in dbase
		self.post_evolve()
		self.db_entry['run_time'] = timer() - b
		self.write_dbase()
	def iteration_output(self):
		output = ( 'gen: ' + `self.gen` + ' ' 
		         + 'max: ' + `self.stats['current']['max']`  + ' ' 
		         + 'dev: ' + `self.p_dev` + ' ' 
		        + 'eval time: ' + `self.step_time` + ' ')
		self._print( output )
			
	def pre_evolve(self): 	pass	
	def post_evolve(self): 	pass	
	def pop_deviation(self):
		#calculate the std deviation across all populations as a percent of mean
		scores = self.pop.scores()
		denom = my_mean(scores)
		if denom == 0.: denom = .0001  # what should I do here?
		return abs(my_std(scores)/denom)
	#dbase stuff
	def init_dbase(self):
		self.db_entry = {}
		self.db_entry['settings'] = self.settings
		t=time.time()
		self.db_entry['raw_time'] = t
		self.db_entry['time'] = time.ctime(t)
		self.db_entry['best_scores'] = [self.stats['current']['max']]
		self.db_entry['stats'] = [copy.deepcopy(self.stats)]
		self.db_entry['step_time'] = [self.step_time]
		self.db_entry['optimization_type'] = string.split(str(self.__class__))[0][1:]
	
	def update_dbase(self):
#		self.db_entry['best_scores'].append(self.pop.best().score())
		self.db_entry['stats'].append(copy.deepcopy(self.stats))
		self.db_entry['step_time'].append(self.step_time)

	def write_dbase(self):	
		"""This does not do file locking on NT - which isn't that big
		   a deal because at the most, two runs are going at a time, and
		   they are unlikely going to write at the same time (but could).
		   On NT, hopefully we're using the gdbm module which does automatic
		   file locking.
		"""
		if(self.settings['dbase'] != ''):
			fname= self.settings['dbase']
			try: 
				if sys.platform == 'win32': pass
				else:
					f = open(fname +'.lock','a')
					fcntl.flock(f.fileno(),fcntl.LOCK_EX)
				try:
					try: db = my_shelve.open(fname,'w')
					except dberror: db = my_shelve.open(fname,'c')	
					keys = db.keys()
					if(len(keys) == 0): self.dbkey = `1`
					else:
						gkeys=[]
						for k in keys:
							try: gkeys.append(string.atoi(k))
							except ValueError: pass
						self.dbkey = `max(gkeys)+1`
					print 'DB NAME: ', self.settings['dbase'], 'KEY: ', self.dbkey
					db[self.dbkey] = self.db_entry 
					db.close()
				except: pass #if an error occured, we still need to unlock the db	
				if sys.platform == 'win32': pass
				else:
					fcntl.flock(f.fileno(),fcntl.LOCK_UN)
					f.close()
			except: 	
				if sys.platform == 'win32': pass
				else:
					f = open('error.lock','a')
					f.write(os.environ['HOST'])
					f.close()

		else:	"no dbase specified"

	def _print(self,val, level = 1):
		if(self.verbose >= level):
			if type(val) == types.StringType: print val
			else:
				pp = pprint.PrettyPrinter(indent=4)
				pp.pprint(val)
	
	
	ALL = -1
class m_galg(galg):
	valid_settings = galg.valid_settings + ['num_pops', 'migrants']
	default_settings = galg.default_settings
	default_settings['pop_size'] = 30; default_settings['num_pops'] = 4
	default_settings['migrants'] = 2
					  
	verbose = 1
	def __init__(self,pop):
		galg.__init__(self,pop)
#		self.GAs = self.GAs + [galg(pop.clone())]
		self.settings = copy.copy(self.default_settings)

	def initialize(self, mode = 'serial'): 
		b = timer()
		#same as galg
		self.test_settings(self.settings)
		self.gen = 0
		sd = self.settings['rand_seed']; alg = self.settings['rand_alg']
		rv.initialize(seed = sd, algorithm = alg)
		self.settings['seed_used'] = rv.initial_seed()
		self._print('initializing... seed = %d' % self.settings['seed_used'])
		self.crossover = self.pop[0].crossover # get the crossover op from the first genome
		self.pop.settings = self.settings
		#end same as galg
		
		#set up my population to hold the best from each sub-pop
		self.pop._size(0) #erase any current member of the pop
		self.pop._size(self.settings['num_pops'])
		self.crossover = self.pop[0].crossover
	
		#extract the galg settings so we don't get a ton of warnings
		#and create the sub ga_s
		sub_ga_settings = {}
		self.GAs = []
		for key in galg.valid_settings:
			sub_ga_settings[key] = self.settings[key]
		for i in range(self.settings['num_pops']): 
			self.GAs.append(galg(self.pop.clone()))
			self.GAs[i].settings = sub_ga_settings.copy()

		self.settings['crossover'] = string.split(str(self.crossover))[0][1:]
		self.settings['selector'] = string.split(str(self.pop.selector))[0][1:]
		self.settings['scaler'] = string.split(str(self.pop.scaler))[0][1:]
		self.settings['genome_type'] = string.split(str(self.pop.model_genome))[0][1:]
		self._print(self.settings)

		if mode[0] == 'p' or mode[0] == 'P':
		    """
			sys.setcheckinterval(1000)
			finished = sync.event()
			bar = sync.barrier(len(self.GAs))
			for ga in self.GAs: 
				thread.start_new_thread(GA_initializer,(bar,finished,ga))
			finished.wait()
			sys.setcheckinterval(10)
			"""
		else:
			for ga in self.GAs: ga.initialize(reseed = 0)								
		cnt = 0		
		for ga in self.GAs: 
			self.pop[cnt] = ga.pop.best()
			cnt = cnt + 1
		self.pop.sort() 			
		self.init_stats()
		self.step_time = timer() - b
		self.init_dbase()

	def init_stats(self):
		#first set up the pops stats, since we don't officially initialize it.
		self.pop.stats = {'current':{},'initial':{},'overall':{}}
		self.pop.stats['selections'] =0; self.pop.stats['crossovers'] =0
		self.pop.stats['mutations'] = 0; self.pop.stats['replacements'] = 0
		self.pop.stats['ind_evals'] = 0
		self.stats = self.pop.stats.copy()
		self.update_stats()
	def update_stats(self):
		"""Gather statistics from the various populations to the mga's population.
		"""
		sum_fields = ['selections','crossovers','mutations','replacements','ind_evals']
		s = []
		for ga in self.GAs:
			for field in sum_fields:
				self.pop.stats[field] = self.pop.stats[field] + ga.stats[field]
			s = s + ga.pop.scores().tolist()			

		self.pop.stats['current']['max'] = self.pop.best().score()
		self.pop.stats['current']['avg'] = my_mean(s)
		self.pop.stats['current']['min'] = min(s)
		if len(s) > 1: self.pop.stats['current']['dev'] = my_std(s)
		else: self.pop.stats['current']['dev'] = 0	
		try: self.pop.stats['overall']['max'] = max(self.pop.stats['overall']['max'],
								    self.pop.stats['current']['max'])
		except KeyError: self.pop.stats['overall']['max'] = self.pop.stats['current']['max']
		try: self.pop.stats['overall']['min'] = min(self.pop.stats['overall']['min'],
								    self.pop.stats['current']['min'])
		except KeyError: self.pop.stats['overall']['min'] = self.pop.stats['current']['min']
		self.pop.stats
		self.pop.stats['pop_evals'] = self.GAs[0].stats['pop_evals']
		self.stats.update(self.pop.stats)


	def step(self,steps=1,mode = 'serial'):
		for st in range(steps):
			b = timer()
			cnt = 0
			#self.pop._size(0) # used if we keep a single pop
			if mode[0] == 'p' or mode[0] == 'P':
			    """
				sys.setcheckinterval(100)
				finished = sync.event()
				bar = sync.barrier(len(self.GAs))
				for ga in self.GAs: 
					thread.start_new_thread(GA_stepper,(bar,finished,ga))
				finished.wait()
				sys.setcheckinterval(10)
				"""
			else:
				for ga in self.GAs: ga.step()								
 
			for ga in self.GAs: 
				#replace the worst member of the local pop 
				self.pop[-1] = ga.pop.best()			
				self.pop.sort()
				#probabaly not the fast approach to things, but... keeps an itelligent pop
				#for ind in ga.pop: self.pop.append(ind)	
	
			self.migrate()
			self.gen = self.gen + 1
			e = timer(); self.step_time = e - b
			self.update_stats()
			self.db_entry['best_scores'].append(self.stats['current']['max'])
			
	def pop_deviation(self):
		"""calculate the std deviation across all populations"""
		all_scores = []
		for ga in self.GAs:
			all_scores = all_scores + ga.pop.scores().tolist()
		if len(all_scores) > 1:
			denom = my_mean(all_scores)
			if denom == 0.: denom = .0001  # what should I do here?
			return abs(my_std(all_scores)/denom)
		return 0

	def evolve(self, mode = 'serial'):
		b = timer()
		self.initialize(mode)
		self.pre_evolve()
		self.p_dev = self.pop_deviation()
		self.iteration_output()
		while (	self.gen < self.settings['gens'] and
				self.settings['p_deviation'] < self.p_dev  ):
			self.step(1,mode)
			self.p_dev = self.pop_deviation()
			self.iteration_output()
		self.update_dbase() #enter status prior to post_evolve in dbase
		self.post_evolve()
		self.db_entry['run_time'] = timer() - b
		self.write_dbase()

	def migrate(self):
		"""Migration moves members from one population to another.  It takes
		   the best N individuals of GAs[0] and puts clones of them into 
		   GAs[1], replacing the worst individuals.  Likewise,
		   GAs[1] best replace GAs[2] worst.  GAs[-1] best are moved
		   to GAs[0].  This 'stepping stone' migration of individuals allows
		   good ideas to move from one population to another, but still
		   allows the individual population to maintain som diversity.
		"""   
		if len(self.GAs) == 1: return
		migrants = self.settings['migrants']
		if migrants > self.settings['pop_size']:
			migrants = self.settings['pop_size']

		movers = []
		for i in range(migrants):
			movers.append(self.GAs[0].pop[i])
		for ga in self.GAs[1:]:	
			for i in range(migrants):
				ga.pop[-i] = movers[i].clone() #replace the worst individual
				movers[i] = ga.pop[i]
		for i in range(migrants):
			self.GAs[0].pop[-i] = movers[i].clone() #replace the worst individual


"""
def GA_stepper(bar,finished,GA):
	t1 = timer()
	GA.step()
	t2 = timer()
	print 'thread ' + `thread.get_ident()` + 'time ' + `t2-t1` + ' sec.'
	bar.enter()
	finished.post()

def GA_initializer(bar,finished,GA):
	t1 = timer()
	GA.initialize(reseed = 0)
	t2 = timer()
	print 'thread ' + `thread.get_ident()` + 'time ' + `t2-t1` + ' sec.'
	bar.enter()
	finished.post()			
"""	
