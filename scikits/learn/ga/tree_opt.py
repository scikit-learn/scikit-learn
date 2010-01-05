"""
    some basic language objects
    if you don't specify a node_type for these, it is usually a good idea to 
    not allow crossovers at these nodes.  That would allow two (maybe different)
    list_range nodes to exchange places.
"""

import tree
import gene


"""
import weakdict
class opt_object(tree.tree_node,weakdict.WeakValue):
	opt_dict = weakdict.WeakDict()
	def __init__(self,node_type, sub_nodes):
		tree.tree_node.__init__(self,sub_nodes,node_type=node_type)
		weakdict.WeakValue.__init__(self)
		opt_object.opt_dict[id(self)] = self
	def _create(self,parent = None):	
		new = tree.tree_node._create(self,parent)
		weakdict.WeakValue.__init__(new)
		self._WeakValue__cid = None #force a reset of the __cid value
		weakdict.WeakValue.__init__(new) #now get a new value for it
		opt_object.opt_dict[id(new)] = new
		return new		
	def __del__(self): 
		tree.tree_node.__del__(self)
		weakdict.WeakValue.__del__(self)
"""		

class opt_object(tree.tree_node):
	def __init__(self,node_type, sub_nodes):
		tree.tree_node.__init__(self,sub_nodes,node_type=node_type)

class float_range(gene.float_gene,opt_object):
	optimize = 1
	def __init__(self,bounds,node_type='float_range', sub_nodes = 0):
		opt_object.__init__(self,node_type,sub_nodes)
		gene.float_gene.__init__(self,bounds[:2])
		tree.tree_node.__init__(self,sub_nodes,node_type=node_type)
		if(len(bounds) == 3): self.default = bounds[2]
		else: self.default = (bounds[0] + bounds[1])/2.
		self._value = self.default
		print self._value
	def clone(self,parent = None):return tree.tree_node.clone(self,parent)			
	"""
	def mutate(self):
		m = gene.float_gene.mutate(self)
		if(m and self.parent): self.parent.recalc(force_parent=1)
		return m
	"""
	def scale(self,sc):
		self.bounds = (self.bounds[0]*sc,self.bounds[1]*sc)
		self.default = self.default*sc
		self._value = self._value*sc
	def defaultize(self):
		self._value = self.default
		for child in self._children: child.defaultize()		
	def create(self,parent):
		new = tree.tree_node.create(self,parent)
		new.initialize() 
		gene.float_gene.initialize(new) 		
		return new
	def __del__(self):
#		gene.float_gene.__del__(self)
		opt_object.__del__(self)	
	def __repr__(self):
		try: 
			val = self.value()
			if( val < .01 or val > 1000): v = "%4.3e" % self.value()
			else: v = "%4.3f" % self.value()				
		except gene.GAError: v = 'not initialized'
		self.label = '%s = %s' % (self.node_type, v)		
		return tree.tree_node.__repr__(self)
import math
class log_float_range(gene.log_float_gene,opt_object):
	optimize = 1
	def __init__(self,bounds,node_type='log_float_range', sub_nodes = 0):
		gene.log_float_gene.__init__(self,bounds[:2])
		opt_object.__init__(self,node_type,sub_nodes)
		if(len(bounds) == 3): self.default = bounds[2]
		else: self.default = (bounds[0] + bounds[1])/2.
		self._value = math.log10(self.default)
	def clone(self,parent = None):return tree.tree_node.clone(self,parent)		
	"""
	def mutate(self):
		m=gene.log_float_gene.mutate(self)
		if(m and self.parent): self.parent.recalc(force_parent=1)
		return m
	"""
	def scale(self,sc):
		self.default = self.default*sc
		sc = log10(sc)
		self.bounds = (self.bounds[0]*sc,self.bounds[1]*sc)
		self._value = self._value*sc
	def defaultize(self):
		self._value = self.default
		for child in self._children: child.defaultize()		
	def create(self,parent):
		new = tree.tree_node.create(self,parent)
		new.initialize() 
		gene.log_float_gene.initialize(new) 		
		return new
	def __del__(self):
#		gene.log_float_gene.__del__(self)
		opt_object.__del__(self)
	def __repr__(self):
		try: 
			val = self.value()
			if( val < .01 or val > 1000): v = "%4.3e" % self.value()
			else: v = "%4.3f" % self.value()				
		except gene.GAError: v = 'not initialized'
		self.label = '%s = %s' % (self.node_type, v)		
		return tree.tree_node.__repr__(self)
		
class list_range(gene.list_gene,opt_object):
	optimize = 1
	def __init__(self,allele_set,node_type='list_range', default=None, sub_nodes = 0):
		gene.list_gene.__init__(self,allele_set)
		opt_object.__init__(self,node_type,sub_nodes)
		gene.list_gene.initialize(self) # prevents trouble in tree generation
		if(default): self.default = default
		else: self.default = allele_set[int(len(allele_set)/2.)] #the center item
		self._value = self.default
	"""
	def mutate(self):
		m=gene.list_gene.mutate(self)
		if(m and self.parent): self.parent.recalc(force_parent=1)
		return m
	"""
	def clone(self,parent = None):return tree.tree_node.clone(self,parent)		
	def scale(self,sc):
		for i in range(len(self.allele_set)): 
			self.allele_set[i] = self.allele_set[i] *sc
		self.default = self.default*sc
		self._value = self._value*sc
	def defaultize(self):
		self._value = self.default
		for child in self._children: child.defaultize()		
	def create(self,parent):
		new = tree.tree_node.create(self,parent)
		new.initialize()
		gene.list_gene.initialize(new) 
		return new
	def __del__(self):
#		gene.list_gene.__del__(self)
		opt_object.__del__(self)
	def __repr__(self):
		self.label = '%s = %s' % (self.node_type, self.value())		
		return tree.tree_node.__repr__(self)

class val(gene.frozen,opt_object):
	optimize = 0
	def __init__(self,val,node_type='val',sub_nodes=0):
		gene.frozen.__init__(self,val)
		opt_object.__init__(self,node_type,sub_nodes)
	def clone(self,parent = None):return tree.tree_node.clone(self,parent)		
	def scale(self,sc): self._value = self._value*sc
	def defaultize(self): pass
	def create(self,parent):
		new = tree.tree_node.create(self,parent)
		new.initialize() 
		return new
	def __del__(self):
#		gene.frozen.__del__(self)
		opt_object.__del__(self)	
	def __repr__(self):
		self.label = '%s = %s' % (self.node_type, self.value())		
		return tree.tree_node.__repr__(self)


"""
These two routines are useful for picking off or replacing the nodes in a 
tree that should be that should be numerically optimized.  They are helpful if
your interested in using a gradient method to optimize some of the paramters
of the array
"""

def pick_numbers(node):
	start = []; lower = []; upper =[];
	for child in node.children():	
		s, l, u = pick_numbers(child)
		start = start + s
		lower = lower + l
		upper = upper + u
	#for now only works with float_genes	
	if hasattr(node,'optimize') and node.optimize == 1: 	
		s = node._value
		l,u = node.bounds
		start = start + [s]
		lower = lower + [l]
		upper = upper + [u]
	else:
		print 'no opt:', node.__class__		
	return start, lower, upper

def put_numbers(node,vals, index = 0):
	for child in node.children():	
		index = put_numbers(child,vals,index)
	if hasattr(node,'optimize') and node.optimize == 1: 	
		s = node._value = vals[index]
		index = index + 1
	return index

"""
Grab the numerical nodes that need to be optimized so that you can directly
manipulate them
"""
def pick_optimize_nodes(node):
	nodes = [];
	for child in node.children():	
		nodes = nodes + pick_optimize_nodes(child)
	#for now only works with float_genes	
	if hasattr(node,'optimize') and node.optimize == 1: 	
		nodes = nodes + [node]
	return nodes