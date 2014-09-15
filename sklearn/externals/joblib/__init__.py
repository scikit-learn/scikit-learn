""" Joblib is a set of tools to provide **lightweight pipelining in
Python**. In particular, joblib offers:

  1. transparent disk-caching of the output values and lazy re-evaluation
     (memoize pattern)

  2. easy simple parallel computing

  3. logging and tracing of the execution

Joblib is optimized to be **fast** and **robust** in particular on large
data and has specific optimizations for `numpy` arrays. It is
**BSD-licensed**.


    ============================== ============================================
    **User documentation**:        http://packages.python.org/joblib

    **Download packages**:         http://pypi.python.org/pypi/joblib#downloads

    **Source code**:               http://github.com/joblib/joblib

    **Report issues**:             http://github.com/joblib/joblib/issues
    ============================== ============================================


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

Joblib strives to address these problems while **leaving your code and
your flow control as unmodified as possible** (no framework, no new
paradigms).

Main features
------------------

1) **Transparent and fast disk-caching of output value:** a memoize or
   make-like functionality for Python functions that works well for
   arbitrary Python objects, including very large numpy arrays. Separate
   persistence and flow-execution logic from domain logic or algorithmic
   code by writing the operations as a set of steps with well-defined
   inputs and  outputs: Python functions. Joblib can save their
   computation to disk and rerun it only if necessary::

      >>> import numpy as np
      >>> from sklearn.externals.joblib import Memory
      >>> mem = Memory(cachedir='/tmp/joblib')
      >>> import numpy as np
      >>> a = np.vander(np.arange(3)).astype(np.float)
      >>> square = mem.cache(np.square)
      >>> b = square(a)                                   # doctest: +ELLIPSIS
      ________________________________________________________________________________
      [Memory] Calling square...
      square(array([[ 0.,  0.,  1.],
             [ 1.,  1.,  1.],
             [ 4.,  2.,  1.]]))
      ___________________________________________________________square - 0...s, 0.0min

      >>> c = square(a)
      >>> # The above call did not trigger an evaluation

2) **Embarrassingly parallel helper:** to make is easy to write readable
   parallel code and debug it quickly::

      >>> from sklearn.externals.joblib import Parallel, delayed
      >>> from math import sqrt
      >>> Parallel(n_jobs=1)(delayed(sqrt)(i**2) for i in range(10))
      [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]


3) **Logging/tracing:** The different functionalities will
   progressively acquire better logging mechanism to help track what
   has been ran, and capture I/O easily. In addition, Joblib will
   provide a few I/O primitives, to easily define define logging and
   display streams, and provide a way of compiling a report.
   We want to be able to quickly inspect what has been run.

4) **Fast compressed Persistence**: a replacement for pickle to work
   efficiently on Python objects containing large data (
   *joblib.dump* & *joblib.load* ).

..
    >>> import shutil ; shutil.rmtree('/tmp/joblib/')

"""

__version__ = '0.8.3'


from .memory import Memory, MemorizedResult
from .logger import PrintTime
from .logger import Logger
from .hashing import hash
from .numpy_pickle import dump
from .numpy_pickle import load
from .parallel import Parallel
from .parallel import delayed
from .parallel import cpu_count
