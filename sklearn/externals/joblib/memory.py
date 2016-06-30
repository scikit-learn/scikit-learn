"""
A context object for caching a function's return value each time it
is called with the same input arguments.

"""

# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org>
# Copyright (c) 2009 Gael Varoquaux
# License: BSD Style, 3 clauses.


from __future__ import with_statement
import os
import shutil
import time
import pydoc
import re
import sys
try:
    import cPickle as pickle
except ImportError:
    import pickle
import functools
import traceback
import warnings
import inspect
import json
import weakref
import io

# Local imports
from . import hashing
from .func_inspect import get_func_code, get_func_name, filter_args
from .func_inspect import format_signature, format_call
from ._memory_helpers import open_py_source
from .logger import Logger, format_time, pformat
from . import numpy_pickle
from .disk import mkdirp, rm_subdirs
from ._compat import _basestring, PY3_OR_LATER

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


def _get_func_fullname(func):
    """Compute the part of part associated with a function.

    See code of_cache_key_to_dir() for details
    """
    modules, funcname = get_func_name(func)
    modules.append(funcname)
    return os.path.join(*modules)


def _cache_key_to_dir(cachedir, func, argument_hash):
    """Compute directory associated with a given cache key.

    func can be a function or a string as returned by _get_func_fullname().
    """
    parts = [cachedir]
    if isinstance(func, _basestring):
        parts.append(func)
    else:
        parts.append(_get_func_fullname(func))

    if argument_hash is not None:
        parts.append(argument_hash)
    return os.path.join(*parts)


def _load_output(output_dir, func_name, timestamp=None, metadata=None,
                 mmap_mode=None, verbose=0):
    """Load output of a computation."""
    if verbose > 1:
        signature = ""
        try:
            if metadata is not None:
                args = ", ".join(['%s=%s' % (name, value)
                                  for name, value
                                  in metadata['input_args'].items()])
                signature = "%s(%s)" % (os.path.basename(func_name),
                                             args)
            else:
                signature = os.path.basename(func_name)
        except KeyError:
            pass

        if timestamp is not None:
            t = "% 16s" % format_time(time.time() - timestamp)
        else:
            t = ""

        if verbose < 10:
            print('[Memory]%s: Loading %s...' % (t, str(signature)))
        else:
            print('[Memory]%s: Loading %s from %s' % (
                    t, str(signature), output_dir))

    filename = os.path.join(output_dir, 'output.pkl')
    if not os.path.isfile(filename):
        raise KeyError(
            "Non-existing cache value (may have been cleared).\n"
            "File %s does not exist" % filename)
    return numpy_pickle.load(filename, mmap_mode=mmap_mode)


# An in-memory store to avoid looking at the disk-based function
# source code to check if a function definition has changed
_FUNCTION_HASHES = weakref.WeakKeyDictionary()


###############################################################################
# class `MemorizedResult`
###############################################################################
class MemorizedResult(Logger):
    """Object representing a cached value.

    Attributes
    ----------
    cachedir: string
        path to root of joblib cache

    func: function or string
        function whose output is cached. The string case is intended only for
        instanciation based on the output of repr() on another instance.
        (namely eval(repr(memorized_instance)) works).

    argument_hash: string
        hash of the function arguments

    mmap_mode: {None, 'r+', 'r', 'w+', 'c'}
        The memmapping mode used when loading from cache numpy arrays. See
        numpy.load for the meaning of the different values.

    verbose: int
        verbosity level (0 means no message)

    timestamp, metadata: string
        for internal use only
    """
    def __init__(self, cachedir, func, argument_hash,
                 mmap_mode=None, verbose=0, timestamp=None, metadata=None):
        Logger.__init__(self)
        if isinstance(func, _basestring):
            self.func = func
        else:
            self.func = _get_func_fullname(func)
        self.argument_hash = argument_hash
        self.cachedir = cachedir
        self.mmap_mode = mmap_mode

        self._output_dir = _cache_key_to_dir(cachedir, self.func,
                                             argument_hash)

        if metadata is not None:
            self.metadata = metadata
        else:
            self.metadata = {}
            # No error is relevant here.
            try:
                with open(os.path.join(self._output_dir, 'metadata.json'),
                          'rb') as f:
                    self.metadata = json.load(f)
            except:
                pass

        self.duration = self.metadata.get('duration', None)
        self.verbose = verbose
        self.timestamp = timestamp

    def get(self):
        """Read value from cache and return it."""
        return _load_output(self._output_dir, _get_func_fullname(self.func),
                            timestamp=self.timestamp,
                            metadata=self.metadata, mmap_mode=self.mmap_mode,
                            verbose=self.verbose)

    def clear(self):
        """Clear value from cache"""
        shutil.rmtree(self._output_dir, ignore_errors=True)

    def __repr__(self):
        return ('{class_name}(cachedir="{cachedir}", func="{func}", '
                'argument_hash="{argument_hash}")'.format(
                    class_name=self.__class__.__name__,
                    cachedir=self.cachedir,
                    func=self.func,
                    argument_hash=self.argument_hash
                    ))

    def __reduce__(self):
        return (self.__class__, (self.cachedir, self.func, self.argument_hash),
                {'mmap_mode': self.mmap_mode})


class NotMemorizedResult(object):
    """Class representing an arbitrary value.

    This class is a replacement for MemorizedResult when there is no cache.
    """
    __slots__ = ('value', 'valid')

    def __init__(self, value):
        self.value = value
        self.valid = True

    def get(self):
        if self.valid:
            return self.value
        else:
            raise KeyError("No value stored.")

    def clear(self):
        self.valid = False
        self.value = None

    def __repr__(self):
        if self.valid:
            return '{class_name}({value})'.format(
                class_name=self.__class__.__name__,
                value=pformat(self.value)
                )
        else:
            return self.__class__.__name__ + ' with no value'

    # __getstate__ and __setstate__ are required because of __slots__
    def __getstate__(self):
        return {"valid": self.valid, "value": self.value}

    def __setstate__(self, state):
        self.valid = state["valid"]
        self.value = state["value"]


###############################################################################
# class `NotMemorizedFunc`
###############################################################################
class NotMemorizedFunc(object):
    """No-op object decorating a function.

    This class replaces MemorizedFunc when there is no cache. It provides an
    identical API but does not write anything on disk.

    Attributes
    ----------
    func: callable
        Original undecorated function.
    """
    # Should be a light as possible (for speed)
    def __init__(self, func):
        self.func = func

    def __call__(self, *args, **kwargs):
        return self.func(*args, **kwargs)

    def call_and_shelve(self, *args, **kwargs):
        return NotMemorizedResult(self.func(*args, **kwargs))

    def __reduce__(self):
        return (self.__class__, (self.func,))

    def __repr__(self):
        return '%s(func=%s)' % (
                    self.__class__.__name__,
                    self.func
            )

    def clear(self, warn=True):
        # Argument "warn" is for compatibility with MemorizedFunc.clear
        pass


###############################################################################
# class `MemorizedFunc`
###############################################################################
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
            numpy arrays. See numpy.load for the meaning of the different
            values.

        compress: boolean, or integer
            Whether to zip the stored data on disk. If an integer is
            given, it should be between 1 and 9, and sets the amount
            of compression. Note that compressed arrays cannot be
            read by memmapping.

        verbose: int, optional
            The verbosity flag, controls messages that are issued as
            the function is evaluated.
    """
    #-------------------------------------------------------------------------
    # Public interface
    #-------------------------------------------------------------------------

    def __init__(self, func, cachedir, ignore=None, mmap_mode=None,
                 compress=False, verbose=1, timestamp=None):
        """
            Parameters
            ----------
            func: callable
                The function to decorate
            cachedir: string
                The path of the base directory to use as a data store
            ignore: list or None
                List of variable names to ignore.
            mmap_mode: {None, 'r+', 'r', 'w+', 'c'}, optional
                The memmapping mode used when loading from cache
                numpy arrays. See numpy.load for the meaning of the
                arguments.
            compress : boolean, or integer
                Whether to zip the stored data on disk. If an integer is
                given, it should be between 1 and 9, and sets the amount
                of compression. Note that compressed arrays cannot be
                read by memmapping.
            verbose: int, optional
                Verbosity flag, controls the debug messages that are issued
                as functions are evaluated. The higher, the more verbose
            timestamp: float, optional
                The reference time from which times in tracing messages
                are reported.
        """
        Logger.__init__(self)
        self.mmap_mode = mmap_mode
        self.func = func
        if ignore is None:
            ignore = []
        self.ignore = ignore

        self._verbose = verbose
        self.cachedir = cachedir
        self.compress = compress
        if compress and self.mmap_mode is not None:
            warnings.warn('Compressed results cannot be memmapped',
                          stacklevel=2)
        if timestamp is None:
            timestamp = time.time()
        self.timestamp = timestamp
        mkdirp(self.cachedir)
        try:
            functools.update_wrapper(self, func)
        except:
            " Objects like ufunc don't like that "
        if inspect.isfunction(func):
            doc = pydoc.TextDoc().document(func)
            # Remove blank line
            doc = doc.replace('\n', '\n\n', 1)
            # Strip backspace-overprints for compatibility with autodoc
            doc = re.sub('\x08.', '', doc)
        else:
            # Pydoc does a poor job on other objects
            doc = func.__doc__
        self.__doc__ = 'Memoized version of %s' % doc

    def _cached_call(self, args, kwargs):
        """Call wrapped function and cache result, or read cache if available.

        This function returns the wrapped function output and some metadata.

        Returns
        -------
        output: value or tuple
            what is returned by wrapped function

        argument_hash: string
            hash of function arguments

        metadata: dict
            some metadata about wrapped function call (see _persist_input())
        """
        # Compare the function code with the previous to see if the
        # function code has changed
        output_dir, argument_hash = self._get_output_dir(*args, **kwargs)
        metadata = None
        # FIXME: The statements below should be try/excepted
        if not (self._check_previous_func_code(stacklevel=4) and
                                 os.path.exists(output_dir)):
            if self._verbose > 10:
                _, name = get_func_name(self.func)
                self.warn('Computing func %s, argument hash %s in '
                          'directory %s'
                        % (name, argument_hash, output_dir))
            out, metadata = self.call(*args, **kwargs)
            if self.mmap_mode is not None:
                # Memmap the output at the first call to be consistent with
                # later calls
                out = _load_output(output_dir, _get_func_fullname(self.func),
                                   timestamp=self.timestamp,
                                   mmap_mode=self.mmap_mode,
                                   verbose=self._verbose)
        else:
            try:
                t0 = time.time()
                out = _load_output(output_dir, _get_func_fullname(self.func),
                                   timestamp=self.timestamp,
                                   metadata=metadata, mmap_mode=self.mmap_mode,
                                   verbose=self._verbose)
                if self._verbose > 4:
                    t = time.time() - t0
                    _, name = get_func_name(self.func)
                    msg = '%s cache loaded - %s' % (name, format_time(t))
                    print(max(0, (80 - len(msg))) * '_' + msg)
            except Exception:
                # XXX: Should use an exception logger
                self.warn('Exception while loading results for '
                          '(args=%s, kwargs=%s)\n %s' %
                          (args, kwargs, traceback.format_exc()))

                shutil.rmtree(output_dir, ignore_errors=True)
                out, metadata = self.call(*args, **kwargs)
                argument_hash = None
        return (out, argument_hash, metadata)

    def call_and_shelve(self, *args, **kwargs):
        """Call wrapped function, cache result and return a reference.

        This method returns a reference to the cached result instead of the
        result itself. The reference object is small and pickeable, allowing
        to send or store it easily. Call .get() on reference object to get
        result.

        Returns
        -------
        cached_result: MemorizedResult or NotMemorizedResult
            reference to the value returned by the wrapped function. The
            class "NotMemorizedResult" is used when there is no cache
            activated (e.g. cachedir=None in Memory).
        """
        _, argument_hash, metadata = self._cached_call(args, kwargs)

        return MemorizedResult(self.cachedir, self.func, argument_hash,
            metadata=metadata, verbose=self._verbose - 1,
            timestamp=self.timestamp)

    def __call__(self, *args, **kwargs):
        return self._cached_call(args, kwargs)[0]

    def __reduce__(self):
        """ We don't store the timestamp when pickling, to avoid the hash
            depending from it.
            In addition, when unpickling, we run the __init__
        """
        return (self.__class__, (self.func, self.cachedir, self.ignore,
                self.mmap_mode, self.compress, self._verbose))

    def format_signature(self, *args, **kwargs):
        warnings.warn("MemorizedFunc.format_signature will be removed in a "
                      "future version of joblib.", DeprecationWarning)
        return format_signature(self.func, *args, **kwargs)

    def format_call(self, *args, **kwargs):
        warnings.warn("MemorizedFunc.format_call will be removed in a "
                      "future version of joblib.", DeprecationWarning)
        return format_call(self.func, args, kwargs)

    #-------------------------------------------------------------------------
    # Private interface
    #-------------------------------------------------------------------------

    def _get_argument_hash(self, *args, **kwargs):
        return hashing.hash(filter_args(self.func, self.ignore,
                                         args, kwargs),
                             coerce_mmap=(self.mmap_mode is not None))

    def _get_output_dir(self, *args, **kwargs):
        """ Return the directory in which are persisted the result
            of the function called with the given arguments.
        """
        argument_hash = self._get_argument_hash(*args, **kwargs)
        output_dir = os.path.join(self._get_func_dir(self.func),
                                  argument_hash)
        return output_dir, argument_hash

    get_output_dir = _get_output_dir  # backward compatibility

    def _get_func_dir(self, mkdir=True):
        """ Get the directory corresponding to the cache for the
            function.
        """
        func_dir = _cache_key_to_dir(self.cachedir, self.func, None)
        if mkdir:
            mkdirp(func_dir)
        return func_dir

    def _hash_func(self):
        """Hash a function to key the online cache"""
        func_code_h = hash(getattr(self.func, '__code__', None))
        return id(self.func), hash(self.func), func_code_h

    def _write_func_code(self, filename, func_code, first_line):
        """ Write the function code and the filename to a file.
        """
        # We store the first line because the filename and the function
        # name is not always enough to identify a function: people
        # sometimes have several functions named the same way in a
        # file. This is bad practice, but joblib should be robust to bad
        # practice.
        func_code = u'%s %i\n%s' % (FIRST_LINE_TEXT, first_line, func_code)
        with io.open(filename, 'w', encoding="UTF-8") as out:
            out.write(func_code)
        # Also store in the in-memory store of function hashes
        is_named_callable = False
        if PY3_OR_LATER:
            is_named_callable = (hasattr(self.func, '__name__')
                                 and self.func.__name__ != '<lambda>')
        else:
            is_named_callable = (hasattr(self.func, 'func_name')
                                 and self.func.func_name != '<lambda>')
        if is_named_callable:
            # Don't do this for lambda functions or strange callable
            # objects, as it ends up being too fragile
            func_hash = self._hash_func()
            try:
                _FUNCTION_HASHES[self.func] = func_hash
            except TypeError:
                # Some callable are not hashable
                pass

    def _check_previous_func_code(self, stacklevel=2):
        """
            stacklevel is the depth a which this function is called, to
            issue useful warnings to the user.
        """
        # First check if our function is in the in-memory store.
        # Using the in-memory store not only makes things faster, but it
        # also renders us robust to variations of the files when the
        # in-memory version of the code does not vary
        try:
            if self.func in _FUNCTION_HASHES:
                # We use as an identifier the id of the function and its
                # hash. This is more likely to falsely change than have hash
                # collisions, thus we are on the safe side.
                func_hash = self._hash_func()
                if func_hash == _FUNCTION_HASHES[self.func]:
                    return True
        except TypeError:
            # Some callables are not hashable
            pass

        # Here, we go through some effort to be robust to dynamically
        # changing code and collision. We cannot inspect.getsource
        # because it is not reliable when using IPython's magic "%run".
        func_code, source_file, first_line = get_func_code(self.func)
        func_dir = self._get_func_dir()
        func_code_file = os.path.join(func_dir, 'func_code.py')

        try:
            with io.open(func_code_file, encoding="UTF-8") as infile:
                old_func_code, old_first_line = \
                            extract_first_line(infile.read())
        except IOError:
                self._write_func_code(func_code_file, func_code, first_line)
                return False
        if old_func_code == func_code:
            return True

        # We have differing code, is this because we are referring to
        # different functions, or because the function we are referring to has
        # changed?

        _, func_name = get_func_name(self.func, resolv_alias=False,
                                     win_characters=False)
        if old_first_line == first_line == -1 or func_name == '<lambda>':
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
        if not old_first_line == first_line and source_file is not None:
            possible_collision = False
            if os.path.exists(source_file):
                _, func_name = get_func_name(self.func, resolv_alias=False)
                num_lines = len(func_code.split('\n'))
                with open_py_source(source_file) as f:
                    on_disk_func_code = f.readlines()[
                        old_first_line - 1:old_first_line - 1 + num_lines - 1]
                on_disk_func_code = ''.join(on_disk_func_code)
                possible_collision = (on_disk_func_code.rstrip()
                                      == old_func_code.rstrip())
            else:
                possible_collision = source_file.startswith('<doctest ')
            if possible_collision:
                warnings.warn(JobLibCollisionWarning(
                        'Possible name collisions between functions '
                        "'%s' (%s:%i) and '%s' (%s:%i)" %
                        (func_name, source_file, old_first_line,
                        func_name, source_file, first_line)),
                                stacklevel=stacklevel)

        # The function has changed, wipe the cache directory.
        # XXX: Should be using warnings, and giving stacklevel
        if self._verbose > 10:
            _, func_name = get_func_name(self.func, resolv_alias=False)
            self.warn("Function %s (stored in %s) has changed." %
                        (func_name, func_dir))
        self.clear(warn=True)
        return False

    def clear(self, warn=True):
        """ Empty the function's cache.
        """
        func_dir = self._get_func_dir(mkdir=False)
        if self._verbose > 0 and warn:
            self.warn("Clearing cache %s" % func_dir)
        if os.path.exists(func_dir):
            shutil.rmtree(func_dir, ignore_errors=True)
        mkdirp(func_dir)
        func_code, _, first_line = get_func_code(self.func)
        func_code_file = os.path.join(func_dir, 'func_code.py')
        self._write_func_code(func_code_file, func_code, first_line)

    def call(self, *args, **kwargs):
        """ Force the execution of the function with the given arguments and
            persist the output values.
        """
        start_time = time.time()
        output_dir, _ = self._get_output_dir(*args, **kwargs)
        if self._verbose > 0:
            print(format_call(self.func, args, kwargs))
        output = self.func(*args, **kwargs)
        self._persist_output(output, output_dir)
        duration = time.time() - start_time
        metadata = self._persist_input(output_dir, duration, args, kwargs)

        if self._verbose > 0:
            _, name = get_func_name(self.func)
            msg = '%s - %s' % (name, format_time(duration))
            print(max(0, (80 - len(msg))) * '_' + msg)
        return output, metadata

    # Make public
    def _persist_output(self, output, dir):
        """ Persist the given output tuple in the directory.
        """
        try:
            mkdirp(dir)
            filename = os.path.join(dir, 'output.pkl')
            numpy_pickle.dump(output, filename, compress=self.compress)
            if self._verbose > 10:
                print('Persisting in %s' % dir)
        except OSError:
            " Race condition in the creation of the directory "

    def _persist_input(self, output_dir, duration, args, kwargs,
                       this_duration_limit=0.5):
        """ Save a small summary of the call using json format in the
            output directory.

            output_dir: string
                directory where to write metadata.

            duration: float
                time taken by hashing input arguments, calling the wrapped
                function and persisting its output.

            args, kwargs: list and dict
                input arguments for wrapped function

            this_duration_limit: float
                Max execution time for this function before issuing a warning.
        """
        start_time = time.time()
        argument_dict = filter_args(self.func, self.ignore,
                                    args, kwargs)

        input_repr = dict((k, repr(v)) for k, v in argument_dict.items())
        # This can fail due to race-conditions with multiple
        # concurrent joblibs removing the file or the directory
        metadata = {"duration": duration, "input_args": input_repr}
        try:
            mkdirp(output_dir)
            with open(os.path.join(output_dir, 'metadata.json'), 'w') as f:
                json.dump(metadata, f)
        except:
            pass

        this_duration = time.time() - start_time
        if this_duration > this_duration_limit:
            # This persistence should be fast. It will not be if repr() takes
            # time and its output is large, because json.dump will have to
            # write a large file. This should not be an issue with numpy arrays
            # for which repr() always output a short representation, but can
            # be with complex dictionaries. Fixing the problem should be a
            # matter of replacing repr() above by something smarter.
            warnings.warn("Persisting input arguments took %.2fs to run.\n"
                          "If this happens often in your code, it can cause "
                          "performance problems \n"
                          "(results will be correct in all cases). \n"
                          "The reason for this is probably some large input "
                          "arguments for a wrapped\n"
                          " function (e.g. large strings).\n"
                          "THIS IS A JOBLIB ISSUE. If you can, kindly provide "
                          "the joblib's team with an\n"
                          " example so that they can fix the problem."
                          % this_duration, stacklevel=5)
        return metadata

    def load_output(self, output_dir):
        """ Read the results of a previous calculation from the directory
            it was cached in.
        """
        warnings.warn("MemorizedFunc.load_output is deprecated and will be "
                      "removed in a future version\n"
                      "of joblib. A MemorizedResult provides similar features",
                      DeprecationWarning)
        # No metadata available here.
        return _load_output(output_dir, _get_func_fullname(self.func),
                            timestamp=self.timestamp,
                            mmap_mode=self.mmap_mode, verbose=self._verbose)

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


###############################################################################
# class `Memory`
###############################################################################
class Memory(Logger):
    """ A context object for caching a function's return value each time it
        is called with the same input arguments.

        All values are cached on the filesystem, in a deep directory
        structure.

        see :ref:`memory_reference`
    """
    #-------------------------------------------------------------------------
    # Public interface
    #-------------------------------------------------------------------------

    def __init__(self, cachedir, mmap_mode=None, compress=False, verbose=1):
        """
            Parameters
            ----------
            cachedir: string or None
                The path of the base directory to use as a data store
                or None. If None is given, no caching is done and
                the Memory object is completely transparent.
            mmap_mode: {None, 'r+', 'r', 'w+', 'c'}, optional
                The memmapping mode used when loading from cache
                numpy arrays. See numpy.load for the meaning of the
                arguments.
            compress: boolean, or integer
                Whether to zip the stored data on disk. If an integer is
                given, it should be between 1 and 9, and sets the amount
                of compression. Note that compressed arrays cannot be
                read by memmapping.
            verbose: int, optional
                Verbosity flag, controls the debug messages that are issued
                as functions are evaluated.
        """
        # XXX: Bad explanation of the None value of cachedir
        Logger.__init__(self)
        self._verbose = verbose
        self.mmap_mode = mmap_mode
        self.timestamp = time.time()
        self.compress = compress
        if compress and mmap_mode is not None:
            warnings.warn('Compressed results cannot be memmapped',
                          stacklevel=2)
        if cachedir is None:
            self.cachedir = None
        else:
            self.cachedir = os.path.join(cachedir, 'joblib')
            mkdirp(self.cachedir)

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
            return functools.partial(self.cache, ignore=ignore,
                                     verbose=verbose, mmap_mode=mmap_mode)
        if self.cachedir is None:
            return NotMemorizedFunc(func)
        if verbose is None:
            verbose = self._verbose
        if mmap_mode is False:
            mmap_mode = self.mmap_mode
        if isinstance(func, MemorizedFunc):
            func = func.func
        return MemorizedFunc(func, cachedir=self.cachedir,
                                   mmap_mode=mmap_mode,
                                   ignore=ignore,
                                   compress=self.compress,
                                   verbose=verbose,
                                   timestamp=self.timestamp)

    def clear(self, warn=True):
        """ Erase the complete cache directory.
        """
        if warn:
            self.warn('Flushing completely the cache')
        if self.cachedir is not None:
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

    def __reduce__(self):
        """ We don't store the timestamp when pickling, to avoid the hash
            depending from it.
            In addition, when unpickling, we run the __init__
        """
        # We need to remove 'joblib' from the end of cachedir
        cachedir = self.cachedir[:-7] if self.cachedir is not None else None
        return (self.__class__, (cachedir,
                self.mmap_mode, self.compress, self._verbose))
