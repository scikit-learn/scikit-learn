from ga_util import *
#import regress
AddError = 'AddError'
import pprint
import sys

#if sys.platform == 'win32':
#	import pywin.debugger
#	pdb = pywin.debugger
#else: 
#	import pdb

#pp = ejpprint.PrettyPrinter(indent = 4)
pp = pprint.PrettyPrinter(indent = 4)

ParentError = 'ParentError'

class base_tree_node:
	objects_ever = 0
	objects = 0
	circular = 0
	
	def inc(self):
		base_tree_node.objects = base_tree_node.objects + 1
		base_tree_node.objects_ever = base_tree_node.objects_ever + 1
	def dec(self):
		base_tree_node.objects = base_tree_node.objects - 1
		
	def __init__(self,child_count,node_type='',derive_type='', parent=None):
#		print 'trenode init',type(parent)
		self.inc()
		if parent is None: self.symbol_table = {}
		self.node_type = node_type
		self.derive_type = derive_type
		self.child_count = child_count
		self._children = []
		self.set_parent(parent)
		self.label = node_type
	#so that we don't overwrite the general creation function in a derived class
	#and screw up clone 	
	def _create(self,parent = None):
		"""
		make lean and mean cause its called a ton
		"""	
		self.inc()
		new = shallow_clone(self)	
		new.set_parent(parent)
		#next 2 lines removed to higher level for speed
		#new._children = [] 
		#if parent is None: new.symbol_table = {}
		return new		
	def create(self,parent = None): 
		new = self._create(parent)
		new._children = []
		if parent is None: new.symbol_table = {}
		return new
	def clone(self,parent = None):
		"""
		make lean and mean cause its called a ton
		"""
		new = self._create(parent) 
		new._children = map(lambda x,par=new: 
					x.clone(par),
					self._children)
						
		if parent is None: new.generate_symbol_table()
		return new
	def set_parent(self,parent):
		self.parent = parent
	def get_parent(self): 
		return self.parent
	def generate_symbol_table(self):
		self.symbol_table = {}
		self._generate_symbol_table(self.symbol_table)
	def _generate_symbol_table(self,symbol_table):
		"""
		Return:
		
		symbol_table -- a dictionary with derive_type
		         symbols as its keys.  Each value is a
		         list of all the nodes in this tree with
		         that derive_type.
		""" 
		for child in self._children:
			child._generate_symbol_table(symbol_table)
		if symbol_table.has_key(self.derive_type):
			symbol_table[self.derive_type].append(self)
		else:	symbol_table[self.derive_type] = [self]
	def parent_test(self):
		for child in self._children:
			child.parent_test()
		for child in self._children:
			if not child.get_parent() is self:
				#pdb.set_trace()
				raise ParentError				
	def node_count(self,type=None):
		cnt = 0
		for child in self._children:
			cnt = cnt + child.node_count(type)
		if self.node_type == type or type is None:
			cnt = cnt + 1
		return cnt	
	
	def apply_node(self,func,args):
		for child in self._children:
			child.apply_node(func,args)
		apply(func,(self,)+args)
		
	def add_child(self,node):
		if len(self._children) < self.child_count:
			node.set_parent(self)
			self._children.append(node)
		else:	raise AddError, 'to many children'
	def filled(self):
		return len(self._children) >= self.child_count
	def children(self):
		return self._children	
	
	def leaves(self):
		if self.child_count == 0:
			return 1
		else:
			leaves = 0
			for child in self._children:
				leaves = leaves + child.leaves()
		return leaves
	def depth(self):
		if self.child_count == 0: return 1
		else:
			return max(map(tree_node.depth,self._children)) + 1
	def ancestors(self):
		if self.get_parent(): return self.get_parent().ancestors() + 1
		return 1
	def root(self):
		if self.get_parent(): return self.get_parent().root()
		return self
	# I needed this to output string chromosomes for antennas
	# rethink this later	
	def file_output(self,file): 
		for child in self._children: child.file_output(file)
	def __repr__(self):
			res = '%s %s' % (self.label, self.derive_type)
			if len(self._children):
#				res = res + ':%s' % self._children 
				res =  res + '\n' + '\t'*self.ancestors()
				res = res + '%s' % self._children 
			return res
	def delete_circulars(self):
		#if hasattr(self,'parent'):
		#	if self.parent is None: print 'deleting root ciculars'
		base_tree_node.circular = base_tree_node.circular + 1
		import sys
		self.symbol_table = None	
		for child in self._children:
			if len(child._children): 
				child.delete_circulars()
			else:
				base_tree_node.circular = base_tree_node.circular + 1
			child.parent = None
			if sys.getrefcount(child) > 3:
			#	print 'aaaah!', sys.getrefcount(child)
				pass
			del child	
	
	def __del__(self):
#		print 'tree_node killed:',tree_node.objects
		self.dec()
		#base_tree_node.objects = base_tree_node.objects - 1
	def __cmp__(self,other):
		#like ga_list compare...
		try: return cmp(self.__dict__,other.__dict__)
		except AttributeError: return 1
		"""
		equal = 0
		try:
			equal = (self.node_type == other.node_type and
				   self.derive_type == other.derive_type and
				   self.child_count == other.child_count and
				   self._children == other._children)
		except AttributeError: pass
		return not equal
		"""
	def __setstate__(self,state):
		for key in state.keys():
			setattr(self, key, state[key])
		self.inc()	
"""			
#core dumps on linux
import weakdict
class weak_tree_node(base_tree_node,weakdict.WeakValue):
	def __init__(self,child_count,node_type='',derive_type='', parent=None):
		weakdict.WeakValue.__init__(self)
		base_tree_node.__init__(self,child_count,node_type,derive_type, parent)
	def set_parent(self,parent):
		print 'in set'
		if not hasattr(self,'parent'):
			self.parent = weakdict.WeakDict()
		if parent: self.parent[0] = parent
		elif self.parent.has_key(0): del self.parent[0]
		print 'out set'
	def get_parent(self): 
		print 'in get'
		if self.parent.has_key(0): p =  self.parent[0]
		else: p = None
		print 'out get'
		return p
	def delete_circulars(self):
		pass
"""
"""
import mxProxy
class proxy_tree_node(base_tree_node):
	passobj = 2 #could be anything
	def set_parent(self,parent):
		self.parent = mxProxy.WeakProxy(parent,None,self.passobj)
	def get_parent(self): 
		if self.parent: return self.parent.proxy_object(self.passobj)
		return None
	def delete_circulars(self):
		pass
"""

tree_node = base_tree_node
#tree_node = weak_tree_node
#tree_node = proxy_tree_node

def ref(): 
	print 'current', base_tree_node.objects
	print 'ever', base_tree_node.objects_ever
	print 'circular deletes', base_tree_node.circular
def test_treenode():
	a = tree_node(2,'root')
	a.add_child(tree_node(0,'kid1'))
	a.add_child(tree_node(0,'kid2'))
	print a
	return a
	"""
#	regress.test('tree_node: add child',a,notes)
	print a

	notes = 'tree_node: add child'
	a = tree_node(2,'root')
	a.add_child(tree_node(0,'kid1'))
	a.add_child(tree_node(0,'kid2'))
	temp = a.children()[0]
	a.children()[0] = a.children()[1]
	a.children()[1] = temp
	print a
	print a.children()[1].root()
	regress.test('tree_node: swap kids',a,notes)
	"""