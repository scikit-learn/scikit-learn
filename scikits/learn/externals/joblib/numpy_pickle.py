"""
A pickler to save numpy arrays in separate .npy files.
"""

# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org> 
# Copyright (c) 2009 Gael Varoquaux
# License: BSD Style, 3 clauses.

import pickle
import traceback
import sys, os

if sys.version_info[0] == 3:
    from pickle import _Unpickler as Unpickler
else:
    from pickle import Unpickler

################################################################################
# Utility objects for persistence.

class NDArrayWrapper(object):
    """ An object to be persisted instead of numpy arrays.

        The only thing this object does, is store the filename in wich
        the array has been persisted.
    """
    def __init__(self, filename):
        self.filename = filename


################################################################################
# Pickler classes

class NumpyPickler(pickle.Pickler):
    """ A pickler subclass that extracts ndarrays and saves them in .npy
        files outside of the pickle.
    """

    def __init__(self, filename):
        self._filename = filename
        self._filenames = [filename, ]
        self.file = open(filename, 'wb')
        # Count the number of npy files that we have created:
        self._npy_counter = 0
        pickle.Pickler.__init__(self, self.file,
                                protocol=pickle.HIGHEST_PROTOCOL)
        # delayed import of numpy, to avoid tight coupling
        import numpy as np
        self.np = np

    def save(self, obj):
        """ Subclass the save method, to save ndarray subclasses in npy
            files, rather than pickling them. Off course, this is a 
            total abuse of the Pickler class.
        """
        if isinstance(obj, self.np.ndarray):
            self._npy_counter += 1
            try:
                filename = '%s_%02i.npy' % (self._filename,
                                            self._npy_counter )
                self._filenames.append(filename)
                self.np.save(filename, obj)
                obj = NDArrayWrapper(os.path.basename(filename))
            except:
                self._npy_counter -= 1
                # XXX: We should have a logging mechanism
                print 'Failed to save %s to .npy file:\n%s' % (
                        type(obj),
                        traceback.format_exc())
        pickle.Pickler.save(self, obj)



class NumpyUnpickler(Unpickler):
    """ A subclass of the Unpickler to unpickle our numpy pickles.
    """
    dispatch = Unpickler.dispatch.copy()

    def __init__(self, filename, mmap_mode=None):
        self._filename = filename
        self.mmap_mode = mmap_mode
        self._dirname  = os.path.dirname(filename)
        self.file = open(filename, 'rb')
        Unpickler.__init__(self, self.file)
        import numpy as np
        self.np = np


    def load_build(self):
        """ This method is called to set the state of a knewly created
            object. 
            
            We capture it to replace our place-holder objects,
            NDArrayWrapper, by the array we are interested in. We
            replace directly in the stack of pickler.
        """
        Unpickler.load_build(self)
        if isinstance(self.stack[-1], NDArrayWrapper):
            nd_array_wrapper = self.stack.pop()
            if self.np.__version__ >= '1.3':
                array = self.np.load(os.path.join(self._dirname,
                                                nd_array_wrapper.filename),
                                                mmap_mode=self.mmap_mode)
            else:
                # Numpy does not have mmap_mode before 1.3
                array = self.np.load(os.path.join(self._dirname,
                                                nd_array_wrapper.filename))
            self.stack.append(array)


    # Be careful to register our new method.
    dispatch[pickle.BUILD] = load_build


################################################################################
# Utility functions

def dump(value, filename):
    """ Persist an arbitrary Python object into a filename, with numpy arrays 
        saved as separate .npy files.

        See Also
        --------
        joblib.load : corresponding loader
    """
    try:
        pickler = NumpyPickler(filename)
        pickler.dump(value)
    finally:
        if 'pickler' in locals() and hasattr(pickler, 'file'):
            pickler.file.flush()
            pickler.file.close()
    return pickler._filenames


def load(filename, mmap_mode=None):
    """ Reconstruct a Python object and the numpy arrays it contains from 
        a persisted file.

        This function loads the numpy array files saved separately. If
        the mmap_mode argument is given, it is passed to np.save and
        arrays are loaded as memmaps. As a consequence, the reconstructed
        object might not match the original pickled object.

        See Also
        --------
        joblib.dump : function to save the object
    """
    try:
        unpickler = NumpyUnpickler(filename, mmap_mode=mmap_mode)
        obj = unpickler.load()
    finally:
        if 'unpickler' in locals() and hasattr(unpickler, 'file'):
            unpickler.file.close()
    return obj

