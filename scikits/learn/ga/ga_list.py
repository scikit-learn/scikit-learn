from ga_util import *
import UserList
import copy

class ga_list(UserList.UserList):
	def data_clone(self):
		new = shallow_clone(self)
		new.data = map(lambda x: x.clone(),self.data)
		return new		
	def touch(self): pass
	def __setslice__(self, i, j, list):
		if type(list) == type(self.data): self.data[i:j] = list
		else: self.data[i:j] = list.data
	def __getslice__(self, i, j):
		new = shallow_clone(self)
		new.data = self.data[i:j]
		new.touch()
		return new
	def __add__(self, list):
		new = shallow_clone(self)
		if type(list) == type(self.data): new.data = self.data + list			
		else: new.data = self.data + list.data			
		new.touch()
		return new
	def __radd__(self, list):
		new = shallow_clone(self)
		if type(list) == type(self.data): new.data = list + self.data		
		else:	new.data = list.data + self.data			
		new.touch()
		return new
	def __mul__(self, n):
		new = shallow_clone(self)
		new.data = self.data*n
		new.touch()
		return new
	__rmul__ = __mul__
	def __cmp__(self, other): return cmp(self.__dict__,other.__dict__)
	
class gene_list(ga_list):
	pass
	"""
	def __setitem__(self, i, item):
 		if(hasattr(item,'_value')): self.data[i] = item
 		else: self.data[i]._value = item
	def __setslice__(self, i, j, list):
		if type(list) == type(self.data):
			if(hasattr(list[0],'_value')): 
				self.data[i:j] = list
			else:
				for k in range(i,j): self.data[k]._value = list[i]	
		else:
			self.data[i:j] = list.data
	"""