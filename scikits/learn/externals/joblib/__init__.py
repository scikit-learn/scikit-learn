""" Joblib is a set of tools to provide **lightweight pipelining in
Python**. In particular, joblib offers:

  1. transparent disk-caching of the output values and lazy re-evaluation
     (memoize pattern)

  2. easy simple parallel computing

  3. logging and tracing of the execution

Joblib is optimized to be fast and robust in particular on large,
long-running functions and has specific optimizations for `numpy` arrays.

____

* The latest user documentation for `joblib` can be found on
  http://packages.python.org/joblib/

* The latest packages can be downloaded from
  http://pypi.python.org/pypi/joblib

* Instructions for developpers can be found at:
  http://github.com/joblib/joblib

joblib is **BSD-licensed**.

Vision
--------

Joblib came out of long-running data-analysis Python scripts. The long
term vision is to provide tools for scientists to achieve better
reproducibility when running jobs, without changing the way numerical
code looks like. However, Joblib can also be used to provide a
light-weight make replacement.

The main problems identified are:

 1) **Lazy evaluation:** People need to rerun over and over the same
    script as it is tuned, but end up commenting out steps, or
    uncommenting steps, as they are needed, as they take long to run.

 2) **Persistence:** It is difficult to persist in an efficient way
    arbitrary objects containing large numpy arrays. In addition,
    hand-written persistence to disk does not link easily the file on
    disk to the corresponding Python object it was persists from in the
    script. This leads to people not a having a hard time resuming the
    job, eg after a crash and persistence getting in the way of work.

The approach taken by Joblib to address these problems is not to build a
heavy framework and coerce user into using it (e.g. with an explicit
pipeline). It strives to leave your code and your flow control as
unmodified as possible.

Current features
------------------

1) **Transparent and fast disk-caching of output value:** a make-like
   functionality for Python functions that works well with large numpy
   arrays. The goal is to separate operations in a set of steps with 
   well-defined inputs and outputs, that are saved and reran only if 
   necessary, by using standard Python functions::

      >>> from joblib import Memory
      >>> mem = Memory(cachedir='/tmp/joblib')
      >>> import numpy as np
      >>> a = np.vander(np.arange(3))
      >>> square = mem.cache(np.square)
      >>> b = square(a)
      ________________________________________________________________________________
      [Memory] Calling square...
      square(array([[0, 0, 1],
             [1, 1, 1],
             [4, 2, 1]]))
      ___________________________________________________________square - 0.0s, 0.0min

      >>> c = square(a)
      >>> # The above call did not trigger an evaluation

2) **Embarrassingly parallel helper:** to make is easy to write readable 
   parallel code and debug it quickly:

      >>> from joblib import Parallel, delayed
      >>> from math import sqrt
      >>> Parallel(n_jobs=1)(delayed(sqrt)(i**2) for i in range(10))
      [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]


3) **Logging/tracing:** The different functionalities will
   progressively acquire better logging mechanism to help track what
   has been ran, and capture I/O easily. In addition, Joblib will
   provide a few I/O primitives, to easily define define logging and
   display streams, and provide a way of compiling a report. 
   We want to be able to quickly inspect what has been run.

Contributing
-------------

The code is `hosted <http://github.com/GaelVaroquaux/joblib>`_ on github.
It is easy to clone the project and experiment with making your own
modifications. If you need extra features, don't hesitate to contribute
them.

.. 
    >>> import shutil ; shutil.rmtree('/tmp/joblib/')

"""

__version__ = '0.4.3'


from .memory import Memory
from .logger import PrintTime, Logger
from .hashing import hash
from .numpy_pickle import dump, load
from .parallel import Parallel, delayed

