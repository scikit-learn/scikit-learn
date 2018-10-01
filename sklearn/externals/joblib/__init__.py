"""Joblib is a set of tools to provide **lightweight pipelining in
Python**. In particular:

1. transparent disk-caching of functions and lazy re-evaluation
   (memoize pattern)

2. easy simple parallel computing

Joblib is optimized to be **fast** and **robust** in particular on large
data and has specific optimizations for `numpy` arrays. It is
**BSD-licensed**.


    ==================== ===============================================
    **Documentation:**       https://joblib.readthedocs.io

    **Download:**            http://pypi.python.org/pypi/joblib#downloads

    **Source code:**         http://github.com/joblib/joblib

    **Report issues:**       http://github.com/joblib/joblib/issues
    ==================== ===============================================


Vision
--------

The vision is to provide tools to easily achieve better performance and
reproducibility when working with long running jobs.

 *  **Avoid computing twice the same thing**: code is rerun over an
    over, for instance when prototyping computational-heavy jobs (as in
    scientific development), but hand-crafted solution to alleviate this
    issue is error-prone and often leads to unreproducible results

 *  **Persist to disk transparently**: persisting in an efficient way
    arbitrary objects containing large data is hard. Using
    joblib's caching mechanism avoids hand-written persistence and
    implicitly links the file on disk to the execution context of
    the original Python object. As a result, joblib's persistence is
    good for resuming an application status or computational job, eg
    after a crash.

Joblib addresses these problems while **leaving your code and your flow
control as unmodified as possible** (no framework, no new paradigms).

Main features
------------------

1) **Transparent and fast disk-caching of output value:** a memoize or
   make-like functionality for Python functions that works well for
   arbitrary Python objects, including very large numpy arrays. Separate
   persistence and flow-execution logic from domain logic or algorithmic
   code by writing the operations as a set of steps with well-defined
   inputs and  outputs: Python functions. Joblib can save their
   computation to disk and rerun it only if necessary::

      >>> from sklearn.externals.joblib import Memory
      >>> cachedir = 'your_cache_dir_goes_here'
      >>> mem = Memory(cachedir)
      >>> import numpy as np
      >>> a = np.vander(np.arange(3)).astype(np.float)
      >>> square = mem.cache(np.square)
      >>> b = square(a)                                   # doctest: +ELLIPSIS
      ________________________________________________________________________________
      [Memory] Calling square...
      square(array([[0., 0., 1.],
             [1., 1., 1.],
             [4., 2., 1.]]))
      ___________________________________________________________square - 0...s, 0.0min

      >>> c = square(a)
      >>> # The above call did not trigger an evaluation

2) **Embarrassingly parallel helper:** to make it easy to write readable
   parallel code and debug it quickly::

      >>> from sklearn.externals.joblib import Parallel, delayed
      >>> from math import sqrt
      >>> Parallel(n_jobs=1)(delayed(sqrt)(i**2) for i in range(10))
      [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]


3) **Fast compressed Persistence**: a replacement for pickle to work
   efficiently on Python objects containing large data (
   *joblib.dump* & *joblib.load* ).

..
    >>> import shutil ; shutil.rmtree(cachedir)

"""

# PEP0440 compatible formatted version, see:
# https://www.python.org/dev/peps/pep-0440/
#
# Generic release markers:
# X.Y
# X.Y.Z # For bugfix releases
#
# Admissible pre-release markers:
# X.YaN # Alpha release
# X.YbN # Beta release
# X.YrcN # Release Candidate
# X.Y # Final release
#
# Dev branch marker is: 'X.Y.dev' or 'X.Y.devN' where N is an integer.
# 'X.Y.dev0' is the canonical version of 'X.Y.dev'
#
__version__ = '0.12.5'


from .memory import Memory, MemorizedResult, register_store_backend
from .logger import PrintTime
from .logger import Logger
from .hashing import hash
from .numpy_pickle import dump
from .numpy_pickle import load
from .compressor import register_compressor
from .parallel import Parallel
from .parallel import delayed
from .parallel import cpu_count
from .parallel import register_parallel_backend
from .parallel import parallel_backend
from .parallel import effective_n_jobs


__all__ = ['Memory', 'MemorizedResult', 'PrintTime', 'Logger', 'hash', 'dump',
           'load', 'Parallel', 'delayed', 'cpu_count', 'effective_n_jobs',
           'register_parallel_backend', 'parallel_backend',
           'register_store_backend', 'register_compressor']
