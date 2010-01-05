"""
Genes are the most basic building block in this genetic algorithm library.
A gene represents a particular trait of an individual solution.  Mutliple genes are
combined together to form a genome.  The entire genome represents a solution
to the problem being solved.

This module contains the base gene class from which all genes are derived.
It also has two gene classes, list_gene and float_gene, that are commonly used.
If you need to create your own gene class, derive it from the class gene.

history:
12/31/98 documentation added ej
"""

from ga_util import *
import scipy.stats as rv
from random import random
import copy

class gene:
    """
    Genes are the most basic building block in this genetic algorithm library.
    A gene represents a particular trait of an individual solution.  The gene class
    contains the current value of a gene object.  It also knows all the possible
    values that the gene can legally have (the allele set).  The gene itself does not 
    know how to initialize itself to a value or mutate to another value.  Instead, it is 
    closely coupled with initializer and mutator classes that perform these duties.
    By changing the initialization and mutation classes, the behavior of a gene can be 
    altered without actually having to create a new gene class.
    
    This class is never actually instantiated.  Instead it is a base class from which 
    specific gene classes are derived.  The duties of the derived class are 
    pretty limited because this class defines almost all the necessary methods.  The 
    derived class must define the allele set for the gene.  It usually will specify 
    an __init__() method, and possibly some other methods specific to the derived
    gene.  It is also necessary to specifiy the initializer and mutator classes 
    to operate on the new gene class
    
    There are several attributes for the class gene:
    
    mutation_rate -- Every gene has a possibility of mutating each generation.
                     This is a value between 0 and 1 that specifies, in percent, 
                     how often the gene mutates.  mutation_rate can either be
                     set for the entire gene class or on a gene by gene basis.
                     Use the set_mutation() function to change the rate
    mr_bounds -- Sometimes it is useful for the mutation rate of the gene to
                 adapt over the course of evolution.  This 2 element tuple specifies
                 the upper and lower bounds that the mutation rate can have.
                 note:  I haven't really messed with this much and think there
                        is probably a better approach
    mutator -- a mutator object (instantiation of a mutator class) that is used 
               to mutate the gene
    initializer -- an intializer object used to initialize the gene
    
    Take a look at float_gene and list_gene to see examples of derived gene classes.
    """
    mr_bounds = (0,.1)
    mutation_rate = .03
    mutator = None
    initializer = None
    is_gene = 1
    def clone(self): 
        """Makes a shallow copy of the object.  override if you need more specialized behavior
        """
        return shallow_clone(self)
    def replicate(self,cnt): 
        """Returns a list with cnt copies of this object in it
        """
        return map(lambda x: x.clone(),[self]*cnt)
    def initialize(self):
        """Calls the initializer objects evaluate() function to initialize the gene 
        """
        self._value = self.initializer.evaluate(self)
        return self.value()
    def set_mutation(self,mrate):
        """
        Set the mutation rate of the gene.
        
            Arguments:

              mrate -- can be one of the following:
              
                * a number between 0 and 1 - sets the mutation rate of the gene to a specific value.
                * "gene" - use the mutation rate set in the class definition for this gene.
                * "adapt" - the mutation rate for the gene is chosen randomly from the range mr_bounds
        """
        if(mrate=='gene'): 
            try: del self.mutation_rate #remove local mrates and use gene classes mrate
            except AttributeError: pass
        elif(mrate=='adapt'): 
            self.mutation_rate = rv.uniform(self.mr_bounds[0],self.mr_bounds[1])[0]
        else: 
            self.__class__.mutation_rate = mrate

    def mutate(self):
        """
            Returns 1 if gene mutated, 0 otherwise.
                
            Calls the **mutator.evaluate()** function to mutate the gene 
            mutation_rate of the time. Otherwise, it does nothing.      
        """ 
        #inlined 'flip_coin' for speed
        if random() < self.mutation_rate: 
            self._value = self.mutator.evaluate(self)
            return 1
        return 0
    def value(self):
        """Return the current value of the gene. """ 
        try: 
            return self._value
        except AttributeError: 
            raise GAError, 'gene not initialized'
    def set_value(self,x):
        """ Set the value of a gene. NO CHECKING!!!
            Don't assign an incompatible value.
        """ 
        self._value = x
    def __repr__(self): 
        try: return `self.value()`
        except GAError: return 'gene not initialized'
    def __add__(self, other):
        try: return self.value() + other.value()
        except AttributeError: return self.value() + other
    __radd__ = __add__
    def __mul__(self, other):
        try: return self.value() * other.value()
        except AttributeError: return self.value() * other
    __rmul__ = __mul__
    def __sub__(self, other):
        try: return self.value() - other.value()
        except AttributeError: return self.value() - other
    def __rsub__(self, other):
        try: return other.value() - self.value()
        except AttributeError: return other - self.value()
    def __div__(self, other):
        try: return self.value() / other.value()
        except: return self.value() / other
    def __rdiv__(self, other):
        try: return other.value() / self.value()
        except AttributeError: return other / self.value()
    def __float__(self): return float(self.value())
    def __complex__(self): return float(self.value())   
    def __neg__(self): return -self.value()
    def __cmp__(self, other):
        try: 
            if self.__class__ == other.__class__ and self.__dict__ == other.__dict__: return 0
        except AttributeError: pass
        v1 = self.value()
        try: v2 = other.value()
        except AttributeError: v2 = other
        return cmp(v1,v2)
        
class list_gene_uniform_mutator:
    """ 
    This class randomly chooses a new gene value from the allele set
    in a list_gene.  It is also useful as an initializer for list_gene.
    """
    def evaluate(self,gene): 
        """ return a randomly chosen value from the genes allele set """
        return rv.choice(gene.allele_set)
class list_gene_gaussian_mutator:
    """ 
    This class chooses a new gene value from the allele set
    in a list_gene.  The new value is chosen from a gaussian 
    distributed distance away from the current values location in the 
    allele set list.  The mutated value is never equal to the current
    gene value.  The dev_width is the standard deviation of the gaussian
    distribution as a percentage of the length of the list.
        
    As an example, suppose a list_gene has the allele_set [0,1,2,3,4,5,6,7,8,9].
    There are 10 entries in this list.  If the dev_width is .1 (the default), 
    then there is a 65% chance the new value will either be 1 position away from
    the current value.  If the current value is 4, then the new value will be
    3 or 5 66% of the time, 2 or 6 29% of the time, and so on based on a gaussian
    distribution.
        
    If the newly chosen index falls outside of the range of the list (for example 
    -1), then a new value is chosen until the value falls inside the lists range.
    The index is NOT truncated to the bottom or top index in the range.    
    """
    def __init__(self,dev_width = .1):
        """Arguments:
            dev_width -- a value between 0 and 1 that specifies the standard
            deviation as a percentage of the length of the list_gene's 
            allele set.
        """
        self.dev_width = dev_width
    def evaluate(self,gene):
        """ return a new value from the genes allele set """
        size = len(gene.allele_set)
        if size == 1: return gene.allele_set[0]
        w = self.dev_width * size
        old = gene.index()
        new = -1; f = -1
        while not (0 <= new < size):
            f = rv.norm(old,w)[0] 
            new = round(f)
            if(old == new and f > new): new = new + 1
            if(old == new and f < new): new = new - 1
        return gene.allele_set[int(new)] 
class list_gene_walk_mutator:
    """ 
      This class chooses a new gene value from the allele set
      in a list_gene.  The newly chosen value is +/-1 element
    in the allele_set from the current gene value. 
    This is like a random walk across the allele_set
    """
    def evaluate(self,gene):
        old = gene.index()
        move = rv.choice((-1,1))
        return gene.allele_set[old + move]
    
class list_gene(gene):
    """
    The value of a list_gene is chosen from a list of
    possible values - the allele_set.
    For example, the gene could be used to represent a
    mathematical oeprator.  Here the allele_set might be
    ['+','-','*','/'].  The list could just as easily be
    a list of numbers (ie. standard capacitor values),
    strings, or anything else.
    
    The default mutator is a gaussian mutator and the 
    default initializer randomly chooses a value from the
    allele_set.
    """
    gaussian_mutator = list_gene_gaussian_mutator
    uniform_mutator = list_gene_uniform_mutator
    walk_mutator = list_gene_walk_mutator
    mutator = gaussian_mutator()
    initializer = uniform_mutator()
    def __init__(self, allele_set): self.allele_set = allele_set
    def index(self,*val):
        """set or retreive a specific value from the allele_set"""
        if len(val): self._value = self.allele_set[val[0]]
        return self.allele_set.index(self.value())  

class list2_gene(list_gene):
    """
    this is something like we'll do to add part variance to capacitor
    and resistor values during evaluation
    """
    func = nop
    def value(self): return func(self._value)
    def __repr__(self): return `self._value` #???

class float_gene_uniform_mutator:
    """ randomly choose a value within the float_gene's bounds"""
    def evaluate(self,gene):
        bounds=gene.bounds
        new =rv.uniform(bounds[0], bounds[1]-bounds[0] ).rvs()[0]
        return new

class float_gene_gaussian_mutator:
    """ 
    chooses a new value for a float_gene with gaussian 
    shaped distribution around the current value.  
    
    dev_width -- a value between 0 and 1.  It is the standard
    deviation for the gaussian distribution as a percentage
    of the float_gene's range.  For example:  If the genes bounds
    are (0,10) and dev_width is .1, then the standard deviation
    is 1.
    """

    def __init__(self,dev_width = .1):
        self.dev_width = dev_width
    def evaluate(self,gene):
        dev = (gene.bounds[1]-gene.bounds[0]) * self.dev_width
        new = gene.bounds[1]
#       while not (gene.bounds[0] <= new < gene.bounds[1]):
#           new = rv.norm(gene.value(),dev)[0]
#       new = rv.norm(gene.value(),dev)[0]
        #get the _value explicitly so mutator will work for log_float also
        new = rv.norm(gene._value,dev).rvs()[0]
        if new > gene.bounds[1]: new = gene.bounds[1]
        if new < gene.bounds[0]: new = gene.bounds[0]
        return new

class float_gene(gene):
    """
    A float_gene is a gene that takes on a floating point value
    between some upper and lower bounds.
    
    The default mutator is a gaussian mutator and the 
    default initializer randomly chooses a value from within
    the upper and lower bounds.
    
    bounds -- A 2 element tuple that specifies the lower and upper
    bounds for the gene.
    """
    gaussian_mutator = float_gene_gaussian_mutator
    uniform_mutator = float_gene_uniform_mutator
    mutator = gaussian_mutator()
    initializer = uniform_mutator()
    def __init__(self,bounds):
        if len(bounds) !=2: raise GAError, 'float_gene: init expects a 2 element tuple of the fomr (min,max)'
        self.bounds = bounds
    def set_value(self,x):
        """ Set the value of a gene. 
            Convertthe value to a float first!
        """ 
        self._value = float(x)

from Numeric import *
from scipy_base.fastumath import *
class log_float_gene(float_gene):
    def __init__(self,bounds):
        if len(bounds) !=2: raise GAError, 'float_gene: init expects a 2 element tuple of the fomr (min,max)'
        self.bounds = log10(array(bounds))
    def value(self):
        """Return the current value of the gene. """ 
        try: return 10.**(self._value)
        except AttributeError: raise GAError, 'gene not initialized'
        
class frozen:
    """frozen is a gene that always maintains the same value.
    """
    def __init__(self,val): self._value = val
    def initialize(self): pass
    def set_mutation(self,mrate): pass
    def mutate(self): pass
    def value(self) : return self._value
    def clone(self): return shallow_clone(self)
    def __float__(self): return float(self._value)
    def __repr__(self): return `self._value`
    def __add__(self, other):
        try: return self._value + other.value()
        except AttributeError: return self._value + other
    __radd__ = __add__
    def __mul__(self, other):
        try: return self._value * other.value()
        except AttributeError: return self._value * other
    __rmul__ = __mul__
    def __sub__(self, other):
        try: return self._value - other.value()
        except AttributeError: return self._value - other
    def __rsub__(self, other):
        try: return other.value() - self._value
        except AttributeError: return other - self._value
    def __div__(self, other):
        try: return self._value / other.value()
        except: return self._value / other
    def __rdiv__(self, other):
        try: return other.value() / self._value
        except AttributeError: return other / self._value
    def __float__(self): return float(self._value)
    def __neg__(self): return -self._value
    def __cmp__(self, other):
        try: 
            if self.__class__ == other.__class__ and self.__dict__ == other.__dict__: return 0
        except AttributeError: pass
        v1 = self.value()
        try: v2 = other.value()
        except AttributeError: v2 = other
        return cmp(v1,v2)

# not sure why this has to be fully qualified, but things are failing otherwise.
# import tree       
from scipy.ga.tree import tree_node
class tree_gene(tree_node):
    mr_bounds = (0,.1)
    mutation_rate = .03
    model_properties = {}
    def __init__(self,child_count,node_type='',derive_type='', parent=None):
        tree_node.__init__(self,child_count,node_type,derive_type, parent)
    def initialize(self,propagate = 1):
        if propagate:
            for child in self._children: child.initialize()         
    def defaultize(self):
        for child in self._children: child.defaultize()     
    def set_mutation(self,mrate):
        if(mrate=='gene'): 
            try: del self.mutation_rate #remove local mrates and use gene classes mrate
            except AttributeError: pass
        elif(mrate=='adapt'): 
            self.mutation_rate = rv.uniform(self.mr_bounds[0],self.mr_bounds[1])[0]
        else: 
            self.__class__.mutation_rate = mrate
        for child in self._children: child.set_mutation(mrate)
                
    def mutate(self,propagate = 1):
        mutated = 0
        #if flip_coin(self.mutation_rate): pass # handle tree mutation
        if propagate:
            for child in self._children: 
                #careful with short circuit "or"
                mutated = child.mutate() or mutated
        return mutated

    def value(self):
        """Return the current value of the gene. """ 
        try: return self._value
        except AttributeError: raise GAError, 'gene not initialized'
    def set_value(self,x):
        """ Set the value of a gene. NO CHECKING!!!
            Don't assign an incompatible value.
        """ 
        self._value = x
