# distutils: language = c
# cython: cdivision = True
# cython: boundscheck = False
# cython: wraparound = False
# cython: profile = False

from _util cimport gcv, ascii_table
import _forward

cdef class Record:

    def __richcmp__(self, other, method):
        if method == 2:
            return self._eq(other)
        elif method == 3:
            return not self._eq(other)
        else:
            return NotImplemented
        
    def _eq(self, other):
        return self.__class__ is other.__class__ and self._getstate() == other._getstate()
        
    def __getitem__(Record self, int idx):
        return self.iterations[idx]
    
    def __len__(Record self):
        return len(self.iterations)
    
    cpdef append(Record self, Iteration iteration):
        self.iterations.append(iteration)
    
    cpdef FLOAT_t mse(Record self, INDEX_t iteration):
        return self.iterations[iteration].get_mse()
    
    cpdef FLOAT_t gcv(Record self, INDEX_t iteration):
        cdef Iteration it = self.iterations[iteration]
        cdef FLOAT_t mse = it.mse
        return gcv(mse,it.get_size(),self.num_samples,self.penalty)
    
    cpdef FLOAT_t rsq(Record self, INDEX_t iteration):
        cdef FLOAT_t mse0 = self.sst#gcv(self.sst,1,self.num_samples,self.penalty)
        cdef FLOAT_t mse = self.mse(iteration)#gcv(self.mse(iteration):,self.iterations[iteration].get_size(),self.num_samples,self.penalty)#self.gcv(iteration)
        return 1 - (mse / mse0)
    
    cpdef FLOAT_t grsq(Record self, INDEX_t iteration):
        cdef FLOAT_t gcv0 = gcv(self.sst,1,self.num_samples,self.penalty)
        cdef FLOAT_t gcv_ = self.gcv(iteration)
        return 1 - (gcv_/gcv0)

cdef class PruningPassRecord(Record):
    def __init__(PruningPassRecord self, INDEX_t num_samples, INDEX_t num_variables, FLOAT_t penalty, FLOAT_t sst, INDEX_t size, FLOAT_t mse):
        self.num_samples = num_samples
        self.num_variables = num_variables
        self.penalty = penalty
        self.sst = sst
        self.iterations = [FirstPruningPassIteration(size, mse)]
        
    def __reduce__(PruningPassRecord self):
        return (PruningPassRecord, (1,1,1.0,1.0,1,1.0), self._getstate())
    
    def _getstate(PruningPassRecord self):
        result = {'num_samples': self.num_samples,
                'num_variables': self.num_variables,
                'penalty': self.penalty,
                'sst': self.sst,
                'iterations': self.iterations,
                'selected': self.selected}
        return result
        
    def __setstate__(PruningPassRecord self, dict state):
        self.num_samples = state['num_samples']
        self.num_variables = state['num_variables']
        self.penalty = state['penalty']
        self.sst = state['sst']
        self.iterations = state['iterations']
        self.selected = state['selected']
        
    cpdef set_selected(PruningPassRecord self, INDEX_t selected):
        self.selected = selected
    
    cpdef INDEX_t get_selected(PruningPassRecord self):
        return self.selected
        
    cpdef roll_back(PruningPassRecord self, Basis basis):
        cdef INDEX_t n = len(self.iterations)
        cdef INDEX_t i
        for i in range(n - self.selected - 1):
            basis[self.iterations[n - i - 1].get_pruned()].unprune()
    
    def __str__(PruningPassRecord self):
        result = ''
        result += 'Pruning Pass\n'
        header = 'iter\tbf\tterms\tmse\tgcv\trsq\tgrsq'.split('\t')
        data = []
        for i, iteration in enumerate(self.iterations):
            row = str(i) + '\t' + str(iteration) + '\t%.3f\t%.3f\t%.3f' % (self.gcv(i),self.rsq(i),self.grsq(i))
            data.append(row.split('\t'))
        result += ascii_table(header,data)
        result += '\nSelected iteration: ' +  str(self.selected) + '\n'
        return result
    
cdef class ForwardPassRecord(Record):
    def __init__(ForwardPassRecord self, INDEX_t num_samples, INDEX_t num_variables, FLOAT_t penalty, FLOAT_t sst):
        self.num_samples = num_samples
        self.num_variables = num_variables
        self.penalty = penalty
        self.sst = sst
        self.iterations = [FirstForwardPassIteration(self.sst)]
        
    def __reduce__(ForwardPassRecord self):
        return (ForwardPassRecord, (1,1,1.0,1.0), self._getstate())
        
    def _getstate(ForwardPassRecord self):
        return {'num_samples': self.num_samples,
                'num_variables': self.num_variables,
                'penalty': self.penalty,
                'sst': self.sst,
                'iterations': self.iterations}
    
    def __setstate__(ForwardPassRecord self, dict state):
        self.num_samples = state['num_samples']
        self.num_variables = state['num_variables']
        self.penalty = state['penalty']
        self.sst = state['sst']
        self.iterations = state['iterations']
        
    cpdef set_stopping_condition(ForwardPassRecord self, int stopping_condition):
        self.stopping_condition = stopping_condition
    
    def __str__(ForwardPassRecord self):
        header = ['iter','parent','var','knot','mse','terms','gcv','rsq','grsq']
        data = []
        for i, iteration in enumerate(self.iterations):
            data.append([str(i)] + str(iteration).split('\t') + ('%.3f\t%.3f\t%.3f' % (self.gcv(i),self.rsq(i),self.grsq(i))).split('\t'))
        result = ''
        result += 'Forward Pass\n'
        result += ascii_table(header, data)
        result += '\nStopping Condition %d: %s\n' % (self.stopping_condition, _forward.stopping_conditions[self.stopping_condition])
        return result

cdef class Iteration:
    
    def __richcmp__(self, other, method):
        if method == 2:
            return self._eq(other)
        elif method == 3:
            return not self._eq(other)
        else:
            return NotImplemented
        
    def _eq(self, other):
        return self.__class__ is other.__class__ and self._getstate() == other._getstate()
    
    cpdef FLOAT_t get_mse(Iteration self):
        return self.mse
    
    cpdef INDEX_t get_size(Iteration self):
        return self.size

cdef class PruningPassIteration(Iteration):
    def __init__(PruningPassIteration self, INDEX_t pruned, INDEX_t size, FLOAT_t mse):
        self.pruned = pruned
        self.size = size
        self.mse = mse
    
    def __reduce__(PruningPassIteration self):
        return (PruningPassIteration, (1,1,1.0), self._getstate())
        
    def _getstate(PruningPassIteration self):
        return {'pruned': self.pruned,
                'size': self.size,
                'mse': self.mse}
    
    def __setstate__(PruningPassIteration self, dict state):
        self.pruned = state['pruned']
        self.size = state['size']
        self.mse = state['mse']
    
    cpdef INDEX_t get_pruned(PruningPassIteration self):
        return self.pruned
        
    def __str__(PruningPassIteration self):
        result = '%s\t%s\t%s' % (str(self.pruned),self.size,'%.2f' % self.mse if self.mse is not None else None)
        return result
    
cdef class FirstPruningPassIteration(PruningPassIteration):
    def __init__(PruningPassIteration self, INDEX_t size, FLOAT_t mse):
        self.size = size
        self.mse = mse
    
    def __reduce__(FirstPruningPassIteration self):
        return (FirstPruningPassIteration, (1,1.0), self._getstate())
    
    def _getstate(FirstPruningPassIteration self):
        return {'size': self.size,
                'mse': self.mse}
    
    def __setstate__(FirstPruningPassIteration self, dict state):
        self.size = state['size']
        self.mse = state['mse']
        
    def __str__(PruningPassIteration self):
        result = '%s\t%s\t%s' % ('-',self.size,'%.2f' % self.mse if self.mse is not None else None)
        return result
    
cdef class ForwardPassIteration(Iteration):
    def __init__(ForwardPassIteration self, INDEX_t parent, INDEX_t variable, int knot, FLOAT_t mse, INDEX_t size):
        self.parent = parent
        self.variable = variable
        self.knot = knot
        self.mse = mse
        self.size = size
        
    def __reduce__(ForwardPassIteration self):
        return (ForwardPassIteration, (1,1,1,1.0,1), self._getstate())
        
    def _getstate(ForwardPassIteration self):
        return {'parent': self.parent,
                'variable': self.variable,
                'knot': self.knot,
                'mse': self.mse,
                'size': self.size}
    
    def __setstate__(ForwardPassIteration self, dict state):
        self.parent = state['parent']
        self.variable = state['variable']
        self.knot = state['knot']
        self.mse = state['mse']
        self.size = state['size']
    
    def __str__(self):
        result = '%d\t%d\t%d\t%4f\t%d' % (self.parent,self.variable,self.knot,self.mse,self.size)
        return result
    
    cpdef set_no_candidates(ForwardPassIteration self, bint value):
        self.no_candidates = value
        
    cpdef no_further_candidates(ForwardPassIteration self):
        return self.no_candidates
    
cdef class FirstForwardPassIteration(ForwardPassIteration):
    def __init__(FirstForwardPassIteration self, FLOAT_t mse):
        self.mse = mse
        
    def __reduce__(FirstForwardPassIteration self):
        return (FirstForwardPassIteration, (1.0,), self._getstate())
    
    def _getstate(FirstForwardPassIteration self):
        return {'mse': self.mse}
    
    def __setstate__(FirstForwardPassIteration self, dict state):
        self.mse = state['mse']
        
    cpdef INDEX_t get_size(FirstForwardPassIteration self):
        return 1
        
    def __str__(self):
        result = '%s\t%s\t%s\t%4f\t%s' % ('-','-','-',self.mse,1)
        return result
    