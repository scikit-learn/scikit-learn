"""
A context object for caching a function's return value each time it
is called with the same input arguments.

"""

# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org> 
# Copyright (c) 2009 Gael Varoquaux
# License: BSD Style, 3 clauses.


import os
import shutil
import sys
import time
import pydoc
try:
    import cPickle as pickle
except ImportError:
    import pickle
import functools
import traceback
import warnings
import inspect
try:
    # json is in the standard library for Python >= 2.6
    import json
except ImportError:
    try:
        import simplejson as json
    except ImportError:
        # Not the end of the world: we'll do without this functionality
        json = None

# Local imports
from .hashing import hash
from .func_inspect import get_func_code, get_func_name, filter_args
from .logger import Logger, format_time
from . import numpy_pickle
from .disk import rm_subdirs

FIRST_LINE_TEXT = "# first line:"

# TODO: The following object should have a data store object as a sub
# object, and the interface to persist and query should be separated in
# the data store.
#
# This would enable creating 'Memory' objects with a different logic for 
# pickling that would simply span a MemorizedFunc with the same
# store (or do we want to copy it to avoid cross-talks?), for instance to
# implement HDF5 pickling. 

# TODO: Same remark for the logger, and probably use the Python logging
# mechanism.

# TODO: Track history as objects are called, to be able to garbage
# collect them.


def extract_first_line(func_code):
    """ Extract the first line information from the function code
        text if available.
    """
    if func_code.startswith(FIRST_LINE_TEXT):
        func_code = func_code.split('\n')
        first_line = int(func_code[0][len(FIRST_LINE_TEXT):])
        func_code = '\n'.join(func_code[1:])
    else:
        first_line = -1
    return func_code, first_line


class JobLibCollisionWarning(UserWarning):
    """ Warn that there might be a collision between names of functions.
    """


################################################################################
# class `Memory`
################################################################################
class MemorizedFunc(Logger):
    """ Callable object decorating a function for caching its return value 
        each time it is called.
    
        All values are cached on the filesystem, in a deep directory
        structure. Methods are provided to inspect the cache or clean it.

        Attributes
        ----------
        func: callable
            The original, undecorated, function.
        cachedir: string
            Path to the base cache directory of the memory context.
        ignore: list or None
            List of variable names to ignore when choosing whether to
            recompute.
        mmap_mode: {None, 'r+', 'r', 'w+', 'c'}
            The memmapping mode used when loading from cache
            numpy arrays. See numpy.load for the meaning of the
            arguments. Only used if save_npy was true when the
            cache was created.
        verbose: int, optional
            The verbosity flag, controls messages that are issued as 
            the function is revaluated.
    """
    #-------------------------------------------------------------------------
    # Public interface
    #-------------------------------------------------------------------------
   
    def __init__(self, func, cachedir, ignore=None, save_npy=True, 
                             mmap_mode=None, verbose=1, timestamp=None):
        """
            Parameters
            ----------
            func: callable
                The function to decorate
            cachedir: string
                The path of the base directory to use as a data store
            ignore: list or None
                List of variable names to ignore.
            save_npy: boolean, optional
                If True, numpy arrays are saved outside of the pickle
                files in the cache, as npy files.
            mmap_mode: {None, 'r+', 'r', 'w+', 'c'}, optional
                The memmapping mode used when loading from cache
                numpy arrays. See numpy.load for the meaning of the
                arguments. Only used if save_npy was true when the
                cache was created.
            verbose: int, optional
                Verbosity flag, controls the debug messages that are issued 
                as functions are revaluated. The higher, the more verbose
            timestamp: float, optional
                The reference time from which times in tracing messages
                are reported.
        """
        Logger.__init__(self)
        self._verbose = verbose
        self.cachedir = cachedir
        self.func = func
        self.save_npy = save_npy
        self.mmap_mode = mmap_mode
        if timestamp is None:
            timestamp = time.time()
        self.timestamp = timestamp
        if ignore is None:
            ignore = []
        self.ignore = ignore
        if not os.path.exists(self.cachedir):
            os.makedirs(self.cachedir)
        try:
            functools.update_wrapper(self, func)
        except:
            " Objects like ufunc don't like that "
        if inspect.isfunction(func):
            doc = pydoc.TextDoc().document(func
                                    ).replace('\n', '\n\n', 1)
        else:
            # Pydoc does a poor job on other objects
            doc = func.__doc__
        self.__doc__ = 'Memoized version of %s' % doc


    def __call__(self, *args, **kwargs):
        # Compare the function code with the previous to see if the
        # function code has changed
        output_dir, _ = self.get_output_dir(*args, **kwargs)
        # FIXME: The statements below should be try/excepted
        if not (self._check_previous_func_code(stacklevel=3) and 
                                 os.path.exists(output_dir)):
            return self.call(*args, **kwargs)
        else:
            try:
                return self.load_output(output_dir)
            except Exception:
                # XXX: Should use an exception logger
                self.warn(
                'Exception while loading results for '
                '(args=%s, kwargs=%s)\n %s' %
                    (args, kwargs, traceback.format_exc())
                    )

                shutil.rmtree(output_dir, ignore_errors=True)
                return self.call(*args, **kwargs)

    #-------------------------------------------------------------------------
    # Private interface
    #-------------------------------------------------------------------------
   
    def _get_func_dir(self, mkdir=True):
        """ Get the directory corresponding to the cache for the
            function.
        """
        module, name = get_func_name(self.func)
        module.append(name)
        func_dir = os.path.join(self.cachedir, *module)
        if mkdir and not os.path.exists(func_dir):
            try:
                os.makedirs(func_dir)
            except OSError:
                """ Dir exists: we have a race condition here, when using 
                    multiprocessing.
                """
                # XXX: Ugly
        return func_dir


    def get_output_dir(self, *args, **kwargs):
        """ Returns the directory in which are persisted the results
            of the function corresponding to the given arguments.

            The results can be loaded using the .load_output method.
        """
        coerce_mmap = (self.mmap_mode is not None)
        argument_hash = hash(filter_args(self.func, self.ignore,
                             *args, **kwargs), 
                             coerce_mmap=coerce_mmap)
        output_dir = os.path.join(self._get_func_dir(self.func),
                                    argument_hash)
        return output_dir, argument_hash
        

    def _write_func_code(self, filename, func_code, first_line):
        """ Write the function code and the filename to a file.
        """
        func_code = '%s %i\n%s' % (FIRST_LINE_TEXT, first_line, func_code)
        file(filename, 'w').write(func_code)


    def _check_previous_func_code(self, stacklevel=2):
        """ 
            stacklevel is the depth a which this function is called, to
            issue useful warnings to the user.
        """
        # Here, we go through some effort to be robust to dynamically
        # changing code and collision. We cannot inspect.getsource
        # because it is not reliable when using IPython's magic "%run".
        func_code, source_file, first_line = get_func_code(self.func)
        func_dir = self._get_func_dir()
        func_code_file = os.path.join(func_dir, 'func_code.py')

        try:
            if not os.path.exists(func_code_file): 
                raise IOError
            old_func_code, old_first_line = \
                            extract_first_line(file(func_code_file).read())
        except IOError:
                self._write_func_code(func_code_file, func_code, first_line)
                return False
        if old_func_code == func_code:
            return True

        # We have differing code, is this because we are refering to
        # differing functions, or because the function we are refering as 
        # changed?

        if old_first_line == first_line == -1:
            _, func_name = get_func_name(self.func, resolv_alias=False,
                                         win_characters=False)
            if not first_line == -1:
                func_description = '%s (%s:%i)' % (func_name, 
                                                source_file, first_line)
            else:
                func_description = func_name
            warnings.warn(JobLibCollisionWarning(
                "Cannot detect name collisions for function '%s'"
                        % func_description), stacklevel=stacklevel)

        # Fetch the code at the old location and compare it. If it is the
        # same than the code store, we have a collision: the code in the
        # file has not changed, but the name we have is pointing to a new
        # code block.
        if (not old_first_line == first_line 
                                    and source_file is not None
                                    and os.path.exists(source_file)):
            _, func_name = get_func_name(self.func, resolv_alias=False)
            num_lines = len(func_code.split('\n'))
            on_disk_func_code = file(source_file).readlines()[
                        old_first_line-1:old_first_line-1+num_lines-1]
            on_disk_func_code = ''.join(on_disk_func_code)
            if on_disk_func_code.rstrip() == old_func_code.rstrip():
                warnings.warn(JobLibCollisionWarning(
                'Possible name collisions between functions '
                "'%s' (%s:%i) and '%s' (%s:%i)" %
                (func_name, source_file, old_first_line, 
                 func_name, source_file, first_line)),
                 stacklevel=stacklevel)

        # The function has changed, wipe the cache directory.
        # XXX: Should be using warnings, and giving stacklevel
        self.clear(warn=True)
        return False


    def clear(self, warn=True):
        """ Empty the function's cache. 
        """
        func_dir = self._get_func_dir(mkdir=False)
        if self._verbose and warn:
            self.warn("Clearing cache %s" % func_dir)
        if os.path.exists(func_dir):
            shutil.rmtree(func_dir, ignore_errors=True)
        try:
            os.makedirs(func_dir)
        except OSError:
            """ Directory exists: it has been created by another process
            in the mean time. """
        func_code, _, first_line = get_func_code(self.func)
        func_code_file = os.path.join(func_dir, 'func_code.py')
        self._write_func_code(func_code_file, func_code, first_line)


    def call(self, *args, **kwargs):
        """ Force the execution of the function with the given arguments and 
            persist the output values.
        """
        start_time = time.time()
        if self._verbose:
            print self.format_call(*args, **kwargs)
        output_dir, argument_hash = self.get_output_dir(*args, **kwargs)
        output = self.func(*args, **kwargs)
        self._persist_output(output, output_dir)
        input_repr = self._persist_input(output_dir, *args, **kwargs)
        duration = time.time() - start_time
        if self._verbose:
            _, name = get_func_name(self.func)
            msg = '%s - %s' % (name, format_time(duration))
            print max(0, (80 - len(msg)))*'_' + msg
        return output


    def format_call(self, *args, **kwds):
        """ Returns a nicely formatted statement displaying the function 
            call with the given arguments.
        """
        path, signature = self.format_signature(self.func, *args,
                            **kwds)
        msg = '%s\n[Memory] Calling %s...\n%s' % (80*'_', path, signature)
        return msg
        # XXX: Not using logging framework
        #self.debug(msg)

    def format_signature(self, func, *args, **kwds):
        # XXX: This should be moved out to a function
        # XXX: Should this use inspect.formatargvalues/formatargspec?
        module, name = get_func_name(func)
        module = [m for m in module if m]
        if module:
            module.append(name)
            module_path = '.'.join(module)
        else:
            module_path = name
        arg_str = list()
        previous_length = 0
        for arg in args:
            arg = self.format(arg, indent=2)
            if len(arg) > 1500:
                arg = '%s...' % arg[:700]
            if previous_length > 80:
                arg = '\n%s' % arg
            previous_length = len(arg)
            arg_str.append(arg)
        arg_str.extend(['%s=%s' % (v, self.format(i)) for v, i in
                                    kwds.iteritems()])
        arg_str = ', '.join(arg_str)

        signature = '%s(%s)' % (name, arg_str)
        return module_path, signature

    # Make make public

    def _persist_output(self, output, dir):
        """ Persist the given output tuple in the directory.
        """
        try:
            if not os.path.exists(dir):
                os.makedirs(dir)
            filename = os.path.join(dir, 'output.pkl')

            if 'numpy' in sys.modules and self.save_npy:
                numpy_pickle.dump(output, filename) 
            else:
                output_file = file(filename, 'w')
                pickle.dump(output, output_file, protocol=2)
                output_file.close()
        except OSError:
            " Race condition in the creation of the directory "


    def _persist_input(self, output_dir, *args, **kwargs):
        """ Save a small summary of the call using json format in the
            output directory.
        """
        argument_dict = filter_args(self.func, self.ignore,
                                    *args, **kwargs)

        input_repr = dict((k, repr(v)) for k, v in argument_dict.iteritems())
        if json is not None:
            # This can fail do to race-conditions with multiple
            # concurrent joblibs removing the file or the directory
            try:
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                json.dump(
                    input_repr,
                    file(os.path.join(output_dir, 'input_args.json'), 'w'),
                    )
            except:
                pass
        return input_repr

    def load_output(self, output_dir):
        """ Read the results of a previous calculation from the directory
            it was cached in.
        """
        if self._verbose > 1:
            t = time.time() - self.timestamp
            print '[Memory]% 16s: Loading %s...' % (
                                    format_time(t),
                                    self.format_signature(self.func)[0]
                                    )
        filename = os.path.join(output_dir, 'output.pkl')
        if self.save_npy:
            return numpy_pickle.load(filename, 
                                     mmap_mode=self.mmap_mode)
        else:
            output_file = file(filename, 'r')
            return pickle.load(output_file)

    # XXX: Need a method to check if results are available.

    #-------------------------------------------------------------------------
    # Private `object` interface
    #-------------------------------------------------------------------------
   
    def __repr__(self):
        return '%s(func=%s, cachedir=%s)' % (
                    self.__class__.__name__,
                    self.func,
                    repr(self.cachedir),
                    )



################################################################################
# class `Memory`
################################################################################
class Memory(Logger):
    """ A context object for caching a function's return value each time it
        is called with the same input arguments.
    
        All values are cached on the filesystem, in a deep directory
        structure.

        see :ref:`memory`
    """
    #-------------------------------------------------------------------------
    # Public interface
    #-------------------------------------------------------------------------
   
    def __init__(self, cachedir, save_npy=True, mmap_mode=None,
                       verbose=1):
        """
            Parameters
            ----------
            cachedir: string or None
                The path of the base directory to use as a data store
                or None. If None is given, no caching is done and
                the Memory object is completely transparent.
            save_npy: boolean, optional
                If True, numpy arrays are saved outside of the pickle
                files in the cache, as npy files.
            mmap_mode: {None, 'r+', 'r', 'w+', 'c'}, optional
                The memmapping mode used when loading from cache
                numpy arrays. See numpy.load for the meaning of the
                arguments. Only used if save_npy was true when the
                cache was created.
            verbose: int, optional
                Verbosity flag, controls the debug messages that are issued 
                as functions are revaluated.
        """
        # XXX: Bad explaination of the None value of cachedir
        Logger.__init__(self)
        self._verbose = verbose
        self.save_npy = save_npy
        self.mmap_mode = mmap_mode
        self.timestamp = time.time()
        if cachedir is None:
            self.cachedir = None
        else:
            self.cachedir = os.path.join(cachedir, 'joblib')
            if not os.path.exists(self.cachedir):
                os.makedirs(self.cachedir)


    def cache(self, func=None, ignore=None, verbose=None,
                        mmap_mode=False):
        """ Decorates the given function func to only compute its return
            value for input arguments not cached on disk.

            Parameters
            ----------
            func: callable, optional
                The function to be decorated
            ignore: list of strings
                A list of arguments name to ignore in the hashing
            verbose: integer, optional
                The verbosity mode of the function. By default that
                of the memory object is used.
            mmap_mode: {None, 'r+', 'r', 'w+', 'c'}, optional
                The memmapping mode used when loading from cache
                numpy arrays. See numpy.load for the meaning of the
                arguments. By default that of the memory object is used.

            Returns
            -------
            decorated_func: MemorizedFunc object
                The returned object is a MemorizedFunc object, that is 
                callable (behaves like a function), but offers extra
                methods for cache lookup and management. See the
                documentation for :class:`joblib.memory.MemorizedFunc`.
        """
        if func is None:
            # Partial application, to be able to specify extra keyword 
            # arguments in decorators
            return functools.partial(self.cache, ignore=ignore)
        if self.cachedir is None:
            return func
        if verbose is None:
            verbose = self._verbose
        if mmap_mode is False:
            mmap_mode = self.mmap_mode
        if isinstance(func, MemorizedFunc):
            func = func.func
        return MemorizedFunc(func, cachedir=self.cachedir,
                                   save_npy=self.save_npy,
                                   mmap_mode=mmap_mode,
                                   ignore=ignore,
                                   verbose=verbose,
                                   timestamp=self.timestamp)


    def clear(self, warn=True):
        """ Erase the complete cache directory.
        """
        if warn:
            self.warn('Flushing completely the cache')
        rm_subdirs(self.cachedir)


    def eval(self, func, *args, **kwargs):
        """ Eval function func with arguments `*args` and `**kwargs`,
            in the context of the memory.

            This method works similarly to the builtin `apply`, except
            that the function is called only if the cache is not
            up to date.

        """
        if self.cachedir is None:
            return func(*args, **kwargs)
        return self.cache(func)(*args, **kwargs)

    #-------------------------------------------------------------------------
    # Private `object` interface
    #-------------------------------------------------------------------------
   
    def __repr__(self):
        return '%s(cachedir=%s)' % (
                    self.__class__.__name__,
                    repr(self.cachedir),
                    )


