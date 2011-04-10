.. _performance-howto:

=========================
How to optimize for speed
=========================

The following gives some practical guidelines to help you write efficient
code for the scikit-learn project.


Python, Cython or C/C++?
========================

In general, the scikit-learn project emphasizes the **readability** of
the source code to make it easy for the project users to dive into the
source code so as to understand how the algorithm behaves on their data
but also for ease of maintanability (by the developers).

When implementing a new algorithm is thus recommended to **start
implementing it in python using numpy and scipy** by taking care of avoiding
looping code using the vectorized idioms of those libraries. In practice
this means trying to **replace any nested for loops by calls to equivalent
numpy array methods**. The goal is to avoid the CPU wasting time in the
python interpreter rather than crunching numbers to fit your statistical
model.

However sometimes an algorithm cannot be expressed efficiently in simple
vectorized numpy code. In this case, the recommended strategy is the
following:

  1. **Profile** the python implementation to find the main bottleneck and isolate
     it in a **dedicated module level function**. This function will be
     reimplemented as a compiled extension module.

  2. If there exists a well maintained BSD or MIT **C/C++** implementation
     of the same algorithm that is not too big, you can write a **cython
     wrapper** for it and include a copy of the source code of the
     library in the scikit-learn source tree: this strategy is used
     for Support Vector Machine and logistic regression (wrappers for
     liblinear and libsvm).

  3. Otherwise, write an optimized version of your python function using
     **cython** directly. This strategy is used for the
     ``scikits.learn.linear_model.ElasticNet`` class for instance.

  4. **Move the python version of the function in the tests** and use it to
     check that the results of the compiled extension are consistent with the
     gold standard, easy to debug python version.

  5. Once the code is optimized (not simple bottleneck spottable by
     profiling), check whether it is possible to have **coarse grained
     parallelism** that is amenable to **multi-processing** by using the
     ``joblib.Parallel`` class.

When using cython, include the generated C source code alongside with
the cython source code. The goal is to make it possible to install the
scikit on any machine with python, numpy, scipy and C/C++ compiler.


.. _profiling-python-code:


Profiling python code
=====================

TODO: sample profiling session with line-prof in ipython

.. _profiling-compiled-extension:


Profiling compiled extensions
=============================

TODO: sample profiling session using YEP


Performance tips for the cython developer
=========================================

TODO: html report, bound checks, division by zero, ...

- http://docs.cython.org/
- http://www.euroscipy.org/file/3696?vid=download
- http://conference.scipy.org/proceedings/SciPy2009/paper_1/
- http://conference.scipy.org/proceedings/SciPy2009/paper_2/


Multi-core parallelism using ``joblib.Parallel``
================================================

TODO: give a simple teaser example here.

Checkout the official joblib documentation:

- http://packages.python.org/joblib/
