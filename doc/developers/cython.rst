.. _cython:

Cython Best Practices, Conventions and Knowledge
================================================

This documents:

* technical aspects of Cython which aren't currently present in Cython's documentation
* tips to develop Cython code in scikit-learn

Tips for developing with Cython in scikit-learn
-----------------------------------------------

Tips to ease development
^^^^^^^^^^^^^^^^^^^^^^^^

* Time spent reading `Cython's documentation <https://cython.readthedocs.io/en/latest/>`_ is not time lost.

* If you intend to use OpenMP: On MacOS, system's distribution of ``clang`` does not implement OpenMP.
  You can install the ``compilers`` package available on ``conda-forge`` which comes with implementations of OpenMP.

* Activating `checks <https://github.com/scikit-learn/scikit-learn/blob/62a017efa047e9581ae7df8bbaa62cf4c0544ee4/sklearn/_build_utils/__init__.py#L68-L87>`_ might help. E.g. for activating boundscheck use:

  .. code-block:: bash

         export SKLEARN_ENABLE_DEBUG_CYTHON_DIRECTIVES=1

* `Start from scratch in a notebook <https://cython.readthedocs.io/en/latest/src/quickstart/build.html#using-the-jupyter-notebook>`_ to understand how to use Cython and to get feedback on your work quickly.

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

* scikit-learn codebase is a bit of a zero of (fused) types (re)definitions, and there's `ongoing work to simplify
  and unify that across the codebase <https://github.com/scikit-learn/scikit-learn/issues/25572>`_. For now, make sure you understand which concrete types are used ultimately.

* It is helpful to use ``gdb`` to debug. In order to do so, one must use a Python interpreter built with debug support
  (debug symbols and proper optimization). To create a new conda environment (which you might need to deactivate and
  reactivate conda after building/installing) with a source-built CPython interpreter:

  .. code-block:: bash

         git clone https://github.com/python/cpython.git
         conda create -n debug-scikit-dev
         conda activate debug-scikit-dev
         cd cpython
         mkdir debug
         cd debug
         ../configure --prefix=$CONDA_PREFIX --with-pydebug
         make EXTRA_CFLAGS='-DPy_DEBUG' -j<num_cores>
         make install

  * You might find this alias handy:

    .. code-block::

         alias cythonX="cython -X language_level=3 -X boundscheck=False -X wraparound=False -X initializedcheck=False -X nonecheck=False -X cdivision=True"

  * which you can use with:

    .. code-block::

         cythonX --annotate source.pyx

Tips for performances
^^^^^^^^^^^^^^^^^^^^^

* Understand the GIL in context for CPython (which problems it solves, what are its limitations) and get a good
  understanding of when Cython will be mapped to C code free of interactions with CPython, when it will not, and when
  it cannot (e.g. presence of interactions with Python objects, which include functions). In this regard,
  `PEP073 <https://peps.python.org/pep-0703/>`_ provides a good overview and context and pathways for removal.

* Make sure you have deactivated `checks <https://github.com/scikit-learn/scikit-learn/blob/62a017efa047e9581ae7df8bbaa62cf4c0544ee4/sklearn/_build_utils/__init__.py#L68-L87>`_.

* Always prefer memoryviews instead over ``cnp.ndarray`` when possible: memoryview are lightweight.

* Avoid memoryview slicing: memoryview slicing might be costly or misleading in some cases and we better not use it IMO,
  even if handling fewer dimensions in some context would be preferable.

* Decorate final classes or methods with ``@final`` (this allows removing virtual tables when needed)

* Inline methods and function when it makes sense

* Make sure your Cython compilation units `use NumPy recent C API <https://github.com/scikit-learn/scikit-learn/blob/62a017efa047e9581ae7df8bbaa62cf4c0544ee4/setup.py#L64-L70>`_.

* In doubt, read the generated C or C++ code if you can: "The fewer C instructions and indirections for a line of
  Cython code, the better" is a good rule of thumb.

Cython Internals
----------------

Concrete method dispatch implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* In most cases, concrete method dispatch is performed **at runtime using virtual tables**
  (i.e. static C-array of pointers to a given types of function, used at runtime to known which
  function to call from an instance). This is costly because there is a double indirection.

* In some cases, concrete method dispatch is performed **statically** (a function is called
  directly by a pointer dereference). This is less costly, especially if there are rooms
  for this function to be inline.

Known Limitations
^^^^^^^^^^^^^^^^^

* No template for types exist in Cython, alternatives are:

  * using fused-types
  * using Tempita, a small preprocessing language to expand Cython source

* Multiple inheritance is not possible for extension types

* No concept similar to interfaces or traits exist

* Currently it is impossible to inherit and override a method which has a fused type in its signature

* No type covariance or type contravariance exist

* Structs cannot be created using constructor syntax with `nogil` `issue and workaround <https://github.com/cython/cython/issues/1642>`_

Performances improvements will come with `Cython 3.0 which has been released soon in beta <https://github.com/cython/cython/issues/4022#issuecomment-1445210880>`_

Basics
^^^^^^

* Understand that ``nogil`` declarations are just hints: when declaring the ``cdef`` functions as nogil,
  means that they can be called without holding the GIL, but it does not release the GIL when entering them.
  You have to do that yourself either by passing ``nogil=True`` to ``cython.parallel.prange`` explicitly,
  or by using an explicit context manager:

.. code-block:: cython

   cdef inline void my_func(self) nogil:

        # Some logic interacting with CPython, e.g. allocating arrays via
        # NumPy.

        with nogil:
            #

        return 0

This item is based on `this comment from Stéfan's Benhel <https://github.com/cython/cython/issues/2798#issuecomment-459971828>`_
