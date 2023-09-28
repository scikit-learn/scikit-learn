.. _cython:

Cython Best Practices, Conventions and Knowledge
================================================

This documents tips to develop Cython code in scikit-learn.

Tips for developing with Cython in scikit-learn
-----------------------------------------------

Tips to ease development
^^^^^^^^^^^^^^^^^^^^^^^^

* Time spent reading `Cython's documentation <https://cython.readthedocs.io/en/latest/>`_ is not time lost.

* If you intend to use OpenMP: On MacOS, system's distribution of ``clang`` does not implement OpenMP.
  You can install the ``compilers`` package available on ``conda-forge`` which comes with an implementation of OpenMP.

* Activating `checks <https://github.com/scikit-learn/scikit-learn/blob/62a017efa047e9581ae7df8bbaa62cf4c0544ee4/sklearn/_build_utils/__init__.py#L68-L87>`_ might help. E.g. for activating boundscheck use:

  .. code-block:: bash

         export SKLEARN_ENABLE_DEBUG_CYTHON_DIRECTIVES=1

* `Start from scratch in a notebook <https://cython.readthedocs.io/en/latest/src/quickstart/build.html#using-the-jupyter-notebook>`_ to understand how to use Cython and to get feedback on your work quickly.
  If you plan to use OpenMP for your implementations in your Jupyter Notebook, do add extra compiler and linkers arguments in the Cython magic.

  .. code-block:: python

         # For GCC and for clang
         %%cython --compile-args=-fopenmp --link-args=-fopenmp
         # For Microsoft's compilers
         %%cython --compile-args=/openmp --link-args=/openmp

* To debug C code (e.g. a segfault), do use ``gdb`` with:

  .. code-block:: bash

         gdb --ex r --args python ./entrypoint_to_bug_reproducer.py

* To have access to some value in place to debug in ``cdef (nogil)`` context, use:

  .. code-block:: cython

         with gil:
             print(state_to_print)

* Note that Cython cannot parse f-strings with ``{var=}`` expressions, e.g.

  .. code-block:: bash

         print(f"{test_val=}")

* scikit-learn codebase has a lot of non-unified (fused) types (re)definitions.
  There currently is `ongoing work to simplify and unify that across the codebase
  <https://github.com/scikit-learn/scikit-learn/issues/25572>`_.
  For now, make sure you understand which concrete types are used ultimately.

* You might find this alias to compile individual Cython extension handy:

    .. code-block::

         # You might want to add this alias to your shell script config.
         alias cythonX="cython -X language_level=3 -X boundscheck=False -X wraparound=False -X initializedcheck=False -X nonecheck=False -X cdivision=True"

         # This generates `source.c` as if you had recompiled scikit-learn entirely.
         cythonX --annotate source.pyx

* Using the ``--annotate`` option with this flag allows generating a HTML report of code annotation.
  This report indicates interactions with the CPython interpreter on a line-by-line basis.
  Interactions with the CPython interpreter must be avoided as much as possible in
  the computationally intensive sections of the algorithms.
  For more information, please refer to `this section of Cython's tutorial <https://cython.readthedocs.io/en/latest/src/tutorial/cython_tutorial.html#primes>`_

    .. code-block::

         # This generates a HTML report (`source.html`) for `source.c`.
         cythonX --annotate source.pyx

Tips for performance
^^^^^^^^^^^^^^^^^^^^

* Understand the GIL in context for CPython (which problems it solves, what are its limitations)
  and get a good understanding of when Cython will be mapped to C code free of interactions with
  CPython, when it will not, and when it cannot (e.g. presence of interactions with Python
  objects, which include functions). In this regard, `PEP073 <https://peps.python.org/pep-0703/>`_
  provides a good overview and context and pathways for removal.

* Make sure you have deactivated `checks <https://github.com/scikit-learn/scikit-learn/blob/62a017efa047e9581ae7df8bbaa62cf4c0544ee4/sklearn/_build_utils/__init__.py#L68-L87>`_.

* Always prefer memoryviews instead over ``cnp.ndarray`` when possible: memoryviews are lightweight.

* Avoid memoryview slicing: memoryview slicing might be costly or misleading in some cases and
  we better not use it, even if handling fewer dimensions in some context would be preferable.

* Decorate final classes or methods with ``@final`` (this allows removing virtual tables when needed)

* Inline methods and function when it makes sense

* Make sure your Cython compilation units `use NumPy recent C API <https://github.com/scikit-learn/scikit-learn/blob/62a017efa047e9581ae7df8bbaa62cf4c0544ee4/setup.py#L64-L70>`_.

* In doubt, read the generated C or C++ code if you can: "The fewer C instructions and indirections
  for a line of Cython code, the better" is a good rule of thumb.

* ``nogil`` declarations are just hints: when declaring the ``cdef`` functions
  as nogil, it means that they can be called without holding the GIL, but it does not release
  the GIL when entering them. You have to do that yourself either by passing ``nogil=True`` to
  ``cython.parallel.prange`` explicitly, or by using an explicit context manager:

    .. code-block:: cython

       cdef inline void my_func(self) nogil:

            # Some logic interacting with CPython, e.g. allocating arrays via NumPy.

            with nogil:
                # The code here is run as is it were written in C.

            return 0

  This item is based on `this comment from Stéfan's Benhel <https://github.com/cython/cython/issues/2798#issuecomment-459971828>`_

* Direct calls to BLAS routines are possible via interfaces defined in ``sklearn.utils._cython_blas``.

Using OpenMP
^^^^^^^^^^^^

Since scikit-learn can be built without OpenMP, it's necessary to protect each
direct call to OpenMP.

The `_openmp_helpers` module, available in
`sklearn/utils/_openmp_helpers.pyx <https://github.com/scikit-learn/scikit-learn/blob/main/sklearn/utils/_openmp_helpers.pyx>`_
provides protected versions of the OpenMP routines. To use OpenMP routines, they
must be ``cimported`` from this module and not from the OpenMP library directly:

.. code-block:: cython

   from sklearn.utils._openmp_helpers cimport omp_get_max_threads
   max_threads = omp_get_max_threads()


The parallel loop, `prange`, is already protected by cython and can be used directly
from `cython.parallel`.
