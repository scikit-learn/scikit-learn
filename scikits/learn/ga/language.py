import types
import  time
from random import *
SymbolError = 'Symbol Error'
DepthError = 'Depth Error'


class language:
	max_depth = 20
	dont_cross = ['ST']  # dont perform crossovers at the start symbol node
	def __init__(self,lang):
		self.lang = lang
		self.dsc = 0
	def generate(self):	
		new = self._gen()
		new.generate_symbol_table()
		return new
	def _gen(self, cur_sym = 'ST', active_node = None, depth = 0):
		#short circuit really large trees
		if(active_node and depth > self.max_depth):
			self.dsc = self.dsc + 1
			active_node.root().delete_circulars()
			raise DepthError
		depth = depth + 1	
		#reset derivation symbol
		if(str(cur_sym) == 'ST' or self.derive_type == ''):
			if(type(cur_sym) == types.StringType):
				self.derive_type = cur_sym
			else:	self.derive_type = cur_sym.node_type
		new_active_node = active_node
		if type(cur_sym) == types.StringType:
			if self.lang.has_key(cur_sym):
				rule = choice(self.lang[cur_sym])
				for sym in rule:
					new_active_node = self._gen(sym,new_active_node,depth)
			else:
				active_node.root().delete_circulars()
				raise SymbolError, "non-terminal symbol not found:%s" % cur_sym
	
		else:
			parent = new_active_node
			new_active_node = cur_sym.create(parent)
			new_active_node.derive_type = self.derive_type
			self.derive_type = ''
			if parent:
				parent.add_child(new_active_node)
			while new_active_node.filled() and new_active_node.get_parent():	
				new_active_node =new_active_node.get_parent()
		return new_active_node

		
