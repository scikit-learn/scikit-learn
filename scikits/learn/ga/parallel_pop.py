from Numeric import *
from scipy_base.fastumath import *

import sys, thread, sync



import remote_exec	

import population



"""

######

#I've got to lean up evaluate and initial in population so that 

#the incorporation of the parallel stuff is smoother.

######

"""

import sys, thread, sync



def array_round(x):

	y = zeros(shape(x))

	for i in range(len(x.flat)): 

		y[i] = int(round(x[i]))

	return y



def divide_list(l,sections):

		Ntot = len(l)

		Nsec = float(sections)

		Neach = Ntot/Nsec

		div_points = array_round(arange(0,Ntot,Neach)).tolist()		

		if div_points[-1] != Ntot: div_points.append(Ntot)

		sub_pops = []

		st = div_points[0]

		for end in div_points[1:]:

			sub_pops.append(l[st:end])

			st = end

		return sub_pops	



class parallel_pop_initializer:

	def evaluate(self,pop,settings = None):

		#only send the individuals out that need evaluation

		if len(pop):

			Nserv = len(pop.server_list)

			groups = divide_list(pop,Nserv)

			sys.setcheckinterval(10)

			finished = sync.event()

			bar = sync.barrier(Nserv)

			print '************',len(groups), len(pop.server_list), len(pop)

			for i in range(len(groups)):

				inputs = {'sub_pop':groups[i],'settings':settings, 'initializer':pop.initializer}

				returns = ('sub_pop',)

				code = 'initializer.evaluate(sub_pop,settings)'

				data_pack = (inputs,returns,code)

				server = pop.server_list[i]

				thread.start_new_thread(remote_thread_init,(bar,finished,server,data_pack))

			finished.wait()

			sys.setcheckinterval(10)

#what is this?			for ind in pop: ind.evaluate(force)



import cPickle

def plen(obj): return len(cPickle.dumps(obj,1))



class parallel_pop_evaluator:

	def evaluate(self,pop,force = 0):

		import tree

		#print '1',tree.ref()

		#only send the individuals out that need evaluation

		if force:	

			_eval_list = pop.data

		else:	

			_eval_list = filter(lambda x: not x.evaluated,pop)

		#print '2',tree.ref()	

		eval_list = pop.clone()

		#print '3',tree.ref()

		eval_list.data = _eval_list

		if len(eval_list):

			Nserv = len(pop.server_list)

			groups = divide_list(eval_list,Nserv)

			#print '4',tree.ref()

			sys.setcheckinterval(10)

			finished = sync.event()

			bar = sync.barrier(Nserv)

			#print "EVAL LENGTH!!!", plen(pop.evaluator)

			gr = groups[0]

			print "GROUP LENGTH!!!", plen(groups[0]), len(gr), 

			#print "IND!!!", plen(gr[0]),plen(gr[0].root)

			#print '4.5',tree.ref()

			for i in range(len(groups)):

				inputs = {'sub_pop':groups[i], 'evaluator':pop.evaluator, 'force':force}

				returns = ('sub_pop',)

				code = 'evaluator.evaluate(sub_pop,force)'

				data_pack = (inputs,returns,code)

				server = pop.server_list[i]

				thread.start_new_thread(remote_thread_eval,(bar,finished,server,data_pack))

			#print '7',tree.ref()	

			finished.wait()

			sys.setcheckinterval(10)

#what is this?			for ind in pop: ind.evaluate(force)

	"""

	def evaluate(self,pop,force = 0):

		#only send the individuals out that need evaluation

		_eval_list = filter(lambda x: not x.evaluated,pop)

		eval_list = pop.clone()

		eval_list.data = _eval_list

		if len(eval_list):

			#finest grain possible

			groups = divide_list(eval_list,len(eval_list))

			finished = sync.event()

			bar = sync.barrier(groups)

			

			sys.setcheckinterval(10)

			Nserv = len(pop.server_list)

			idx = 0

			while idx < len(groups):

				inputs = {'sub_pop':groups[idx], 'evaluator':pop.evaluator}

				returns = ('sub_pop',)

				code = 'evaluator.evaluate(sub_pop)'

				data_pack = (inputs,returns,code)

				server = pop.server_list[i]

				thread.start_new_thread(remote_thread_eval,(bar,finished,server,data_pack))

			#for i in range(len(groups)):

			#	inputs = {'sub_pop':groups[i], 'evaluator':pop.evaluator}

			#	returns = ('sub_pop',)

			#	code = 'evaluator.evaluate(sub_pop)'

			#	data_pack = (inputs,returns,code)

			#	server = pop.server_list[i]

			#	thread.start_new_thread(remote_thread,(bar,finished,server,data_pack))

			finished.wait()

			sys.setcheckinterval(10)

#what is this?			for ind in pop: ind.evaluate(force)

	"""



def remote_thread_init(bar,finished,server,data_pack):

	try:

		remote = remote_exec.remote_exec(server[0],server[1],0,1)

		results = remote.run(data_pack)

		#assign the results from the returned data to the local individuals

		inputs = data_pack[0]

		old = inputs['sub_pop']

		new = results['sub_pop']

		for i in range(len(old)):

			old[i].__dict__.update(new[i].__dict__)

	except IndexError:

		print 'error in %s,%d' %  server

	bar.enter()

	finished.post()



def remote_thread_eval(bar,finished,server,data_pack):

	import tree

	try:

		#print '5',tree.ref()

		remote = remote_exec.remote_exec(server[0],server[1],0,1)

		results = remote.run(data_pack)

		#print '6',tree.ref()

		#assign the results from the returned data to the local individuals

		inputs = data_pack[0]

		old = inputs['sub_pop']

		new = results['sub_pop']

		for gnm in new:

			gnm.root.delete_circulars()

			del gnm.root

		#print '6.25',tree.ref()

		for i in range(len(old)):

			old[i].__dict__.update(new[i].__dict__)



		#print '6.5',tree.ref()

	except IndexError:

		print 'error in %s,%d' %  server

	"""

	import sys

	#r = new[0].root

	#print 'ref count',sys.getrefcount(r)

	#print '6.75',tree.ref()		

	#Huh??? Why do I need to delete the new genomes

	#individually here?  Why aren't they garbage collected?

	indices = range(len(new))

	indices.reverse()

	for i in indices:

		del new[i]

	#print 'ref count',sys.getrefcount(r)

	#print '6.8',tree.ref()	

	#r.delete_circulars()	

	#print 'ref count',sys.getrefcount(r)

	#print '6.9',tree.ref()	

	#del r

	#print '6.95',tree.ref()	

	"""

	bar.enter()

	finished.post()



class ga_parallel_pop(population.population):

	parallel_evaluator = parallel_pop_evaluator()

	parallel_initializer = parallel_pop_initializer()

	def __init__(self,genome,size=1,server_list=None):

		"""Arguments:

		

		   genome -- a genome object.

		   size -- number.  The population size.  The genome will be 

		           replicated size times to fill the population.

		   server_list -- a list of tuple pairs with machine names and

		                  ports listed for the available servers

		                  ex: [(ee.duke.edu,8000),('elsie.ee.duke.edu',8000)]        

		"""

		population.population.__init__(self,genome,size)

		assert(server_list)

		self.server_list = server_list

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

		self.parallel_initializer.evaluate(self,settings)

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

	

	def evaluate(self, force = 0):

		""" call the parallel_evaluator instead of the evaluator directly

		"""

		self.selector.clear()

		self.parallel_evaluator.evaluate(self,force)

		#self.post_evaluate()

		#all of the remaining should be put in post eval...

		self.sort()

		#this is a cluge to get eval count to work correctly

		preval = self.stats['ind_evals']

		for ind in self:  

			self.stats['ind_evals'] = self.stats['ind_evals'] + ind.evals

			ind.evals = 0		

		print 'evals: ', self.stats['ind_evals'] - preval

		self.touch()

		self.evaluated = 1



########################## test stuff ############################

import genome

import gene

import time



import socket			



class objective:

	def __init__(self,wait=.01):

		self.wait = wait

	def evaluate(self,genome): 

		time.sleep(self.wait)

		return sum(genome.array())



def test_pop(server_list,size=100,wait=.01):

	obj = objective(wait)

	the_gene = gene.float_gene((0,2.5))

	genome = genome.list_genome(the_gene.replicate(5))

	genome.evaluator = obj

	pop = ga_parallel_pop(genome,size,server_list)

	print  '########### awaiting evaluation#############'

	pop.initialize()

	print ' evaluation done!'

	print 'best:', pop.best()

	print 'worst',pop.worst()



def gen_pop():

	genome.list_genome.evaluator = objective()

	gene = gene.float_gene((0,2.5))

	genome = genome.list_genome(gene.replicate(5))

	pop = ga_parallel_pop(genome,100,[(host,port),])

	return pop



	import parallel_pop,beowulf,os



def test_pop2(server_list,size=100,wait=.01):

	import hmm_gnm,os

	genome = hmm_gnm.make_genome()

	#pop = ga_parallel_pop(genome,4,server_list)

	global galg

	#genome.target = targets[0]

	pop = ga_parallel_pop(genome,1,server_list)

	galg = hmm_gnm.class_ga(pop)

	galg.settings.update({ 'pop_size':6,'gens':2,'p_mutate':.03,

				    'dbase':os.environ['HOME'] + '/all_lift3', 'p_cross':0.9, 'p_replace':.6,

				    'p_deviation': -.001})

	galg.evolve()



	print  '########### awaiting evaluation#############'

	pop.initialize()

	print ' evaluation done!'

	print 'best:', pop.best()

	print 'worst',pop.worst()

	

import thread

def test():

	host = socket.gethostname()

	port = 8000

	server_list = [(host,port),(host,port+1)]

	for server in server_list:

		host,port = server

		thread.start_new_thread(remote_exec.server,(host,port))

	thread.start_new_thread(test_pop2,(server_list,))		

	

def test2(machines=32,size=100,wait=.01):

	import time

	t1 = time.time()

	#requires that servers are started on beowulf 1 and 2.

	import beowulf

	server_list = beowulf.beowulf.servers[:machines]

	thread.start_new_thread(test_pop,(server_list,size,wait))			

	print 'total time:', time.time()-t1
