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
     :class:`scikits.learn.linear_model.ElasticNet` class for instance.

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

In order to profile python we recommend to write a script that loads
and prepare you data and then use the ipython integrated python profiler
for interactively exploring the relevant part for the code.

Suppose we want to profile the Non Negative Matrix Factorization module
of the scikit. Let us setup a new ipython session and load the digits dataset
and as in the :ref:`example_decomposition_plot_nmf.py` example::

  In [1]: from scikits.learn.decomposition import NMF

  In [2]: from scikits.learn.datasets import load_digits

  In [3]: X = load_digits().data

Let us first have a look at the overall performance profile using the ``%prun``
magic command::

  In [4]: %prun -l nmf.py NMF(n_components=16, tol=1e-2).fit(X)
           14496 function calls in 1.682 CPU seconds

     Ordered by: internal time
     List reduced from 90 to 9 due to restriction <'nmf.py'>

     ncalls  tottime  percall  cumtime  percall filename:lineno(function)
         36    0.609    0.017    1.499    0.042 nmf.py:151(_nls_subproblem)
       1263    0.157    0.000    0.157    0.000 nmf.py:18(_pos)
          1    0.053    0.053    1.681    1.681 nmf.py:352(fit_transform)
        673    0.008    0.000    0.057    0.000 nmf.py:28(norm)
          1    0.006    0.006    0.047    0.047 nmf.py:42(_initialize_nmf)
         36    0.001    0.000    0.010    0.000 nmf.py:36(_sparseness)
         30    0.001    0.000    0.001    0.000 nmf.py:23(_neg)
          1    0.000    0.000    0.000    0.000 nmf.py:337(__init__)
          1    0.000    0.000    1.681    1.681 nmf.py:461(fit)


Note the use of the ``-l nmf.py`` that restricts the output to lines that
contains the "nmf.py" string. This is useful to have a quick look at the hotspot
of the nmf python module it-self ignoring anything else.

The above results show that the ``_nls_subproblem`` function is the hotspot: it
takes almost 70% of the time of the module. In order to better understand the
profile of this specific function, let us install ``line-prof`` and wire it to
ipython::

  $ pip install line-prof
  $ vim ~/.ipython/ipy_user_conf.py

Ensure the following lines are present::

  import IPython.ipapi
  ip = IPython.ipapi.get()

Towards the end of the file, define the ``%lprun`` magic::

  import line_profiler
  ip.expose_magic('lprun', line_profiler.magic_lprun)

Now restart ipython and let us use this new toy::

  In [1]: from scikits.learn.datasets import load_digits

  In [2]: from scikits.learn.decomposition.nmf import _nls_subproblem, NMF

  In [3]: X = load_digits().data

  In [4]: %lprun -f _nls_subproblem NMF(n_components=16, tol=1e-2).fit(X)
  Timer unit: 1e-06 s

  File: scikits/learn/decomposition/nmf.py
  Function: _nls_subproblem at line 137
  Total time: 1.73153 s

  Line #      Hits         Time  Per Hit   % Time  Line Contents
  ==============================================================
     137                                           def _nls_subproblem(V, W, H_init, tol, max_iter):
     138                                               """Non-negative least square solver
     ...
     170                                               """
     171        48         5863    122.1      0.3      if (H_init < 0).any():
     172                                                   raise ValueError("Negative values in H_init passed to NLS solver.")
     173
     174        48          139      2.9      0.0      H = H_init
     175        48       112141   2336.3      5.8      WtV = np.dot(W.T, V)
     176        48        16144    336.3      0.8      WtW = np.dot(W.T, W)
     177
     178                                               # values justified in the paper
     179        48          144      3.0      0.0      alpha = 1
     180        48          113      2.4      0.0      beta = 0.1
     181       638         1880      2.9      0.1      for n_iter in xrange(1, max_iter + 1):
     182       638       195133    305.9     10.2          grad = np.dot(WtW, H) - WtV
     183       638       495761    777.1     25.9          proj_gradient = norm(grad[np.logical_or(grad < 0, H > 0)])
     184       638         2449      3.8      0.1          if proj_gradient < tol:
     185        48          130      2.7      0.0              break
     186
     187      1474         4474      3.0      0.2          for inner_iter in xrange(1, 20):
     188      1474        83833     56.9      4.4              Hn = H - alpha * grad
     189                                                       # Hn = np.where(Hn > 0, Hn, 0)
     190      1474       194239    131.8     10.1              Hn = _pos(Hn)
     191      1474        48858     33.1      2.5              d = Hn - H
     192      1474       150407    102.0      7.8              gradd = np.sum(grad * d)
     193      1474       515390    349.7     26.9              dQd = np.sum(np.dot(WtW, d) * d)
     ...

By looking at the top values of the ``% Time`` column it is really easy to
pin-point the most expensive expressions that would deserve additional care.


.. _profiling-compiled-extension:

Profiling compiled extensions
=============================

TODO: sample profiling session using YEP

- https://github.com/fabianp/yep
- http://fseoane.net/blog/2011/a-profiler-for-python-extensions/


Performance tips for the cython developer
=========================================

TODO: html report, type declarations, bound checks, division by zero checks,
memory alignement, direct blas calls...

- http://docs.cython.org/
- http://www.euroscipy.org/file/3696?vid=download
- http://conference.scipy.org/proceedings/SciPy2009/paper_1/
- http://conference.scipy.org/proceedings/SciPy2009/paper_2/


Multi-core parallelism using ``joblib.Parallel``
================================================

TODO: give a simple teaser example here.

Checkout the official joblib documentation:

- http://packages.python.org/joblib/
