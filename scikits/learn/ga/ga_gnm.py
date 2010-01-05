"""

This module adds gradient optimization capabilities to the 

standard genomes.  Basically this means that any problem set up

for a GA is automatically able to be gradient optimized...



For purity sake, the grad and genome

modules have been left totally separate.  It might have

been just as easy to derive genomes directly from 

grad.grad - and maybe that will happen in the future.



Caveats:

	This has only be set up for list_genomes made up of 

	floating point genes.  The tree_genomes just need to

	be recoded here translating the pick_numbers functions

	from tree_opt.  

	

	genomes of discrete variable genes should be able to work also.

"""



import grad

import genome



class list_genome(genome.list_genome,grad.grad):

	""" So far, grad_min, and grad_max only 

		work for float_genes.

		Test:

		#Test gradient optimization

		>>> import ga_gnm, gene

		>>> g = gene.float_gene((-1,1))

		>>> class simple_genome(ga_gnm.list_genome):

		...		def performance(self):

		...			s = 0

		...			for i in self: s = s+ i

		...			return s							

		>>> a = simple_genome(g.replicate(10))

		>>> a.initialize()

		>>> a.grad_opt(5)

		33

		>>> a

		[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

	"""

	def grad_params(self):

		return self.get_values() # calls list__genome	get_values()

	def set_grad_params(self,x):

		self.set_values(x) # calls list__genome	set_values()

	#do we really need this?	

	def grad_len(self):

		return len(self)

	def grad_min(self):

		gmin = []

		for flt_gene in self: 

			gmin.append(flt_gene.bounds[0])	

		return gmin	

	def grad_max(self):

		gmax = []

		for flt_gene in self: 

			gmax.append(flt_gene.bounds[1])	

		return gmax	

