.. _install_bleeding_edge:

Set up your development environment
-----------------------------------

.. _git_repo:

Fork the scikit-learn repository
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
First, you need to `create an account <https://github.com/join>`_ on
GitHub (if you do not already have one) and fork the `project repository
<https://github.com/scikit-learn/scikit-learn>`__ by clicking on the 'Fork'
button near the top of the page. This creates a copy of the code under your
account on the GitHub user account. For more details on how to fork a
repository see `this guide <https://help.github.com/articles/fork-a-repo/>`_.

The following steps explain how to set up a local clone of your forked git repository
and how to locally install scikit-learn according to your operating system.

Set up a local clone of your fork
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Clone your fork of the scikit-learn repo from your GitHub account to your
local disk:

.. prompt:: bash

  git clone https://github.com/YourLogin/scikit-learn.git  # add --depth 1 if your connection is slow

and change into that directory:

.. prompt:: bash

  cd scikit-learn

.. _upstream:

Next, add the ``upstream`` remote. This saves a reference to the main
scikit-learn repository, which you can use to keep your repository
synchronized with the latest changes (you'll need this later in the :ref:`development_workflow`):

.. prompt:: bash

  git remote add upstream https://github.com/scikit-learn/scikit-learn.git

Check that the `upstream` and `origin` remote aliases are configured correctly
by running:

.. prompt:: bash

  git remote -v

This should display:

.. code-block:: text

  origin    https://github.com/YourLogin/scikit-learn.git (fetch)
  origin    https://github.com/YourLogin/scikit-learn.git (push)
  upstream  https://github.com/scikit-learn/scikit-learn.git (fetch)
  upstream  https://github.com/scikit-learn/scikit-learn.git (push)


Set up a dedicated environment and install dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
..
   TODO Add |PythonMinVersion| to min_dependency_substitutions.rst one day.
   Probably would need to change a bit sklearn/_min_dependencies.py since Python is not really a package ...
.. |PythonMinVersion| replace:: 3.10

Using an isolated environment such as `venv` or `conda` makes it possible to
install a specific version of scikit-learn with pip or conda and its dependencies
independently of any previously installed Python packages. Note that the virtual
environment is optional but strongly recommended, in order to avoid potential
conflicts with other packages.

In addition to the required Python dependencies, you need to have a working C/C++
compiler with OpenMP_ support to build scikit-learn Cython extensions.

.. note::

      If OpenMP is not supported by the compiler, the build will be done with
      OpenMP functionalities disabled. This is not recommended since it will force
      some estimators to run in sequential mode instead of leveraging thread-based
      parallelism. Setting the ``SKLEARN_FAIL_NO_OPENMP`` environment variable
      (before cythonization) will force the build to fail if OpenMP is not
      supported.

In the following, you find the specific install instructions for all supported platforms and package managers.

.. raw:: html

  <style>
    /* Show caption on large screens */
    @media screen and (min-width: 960px) {
      .install-instructions .sd-tab-set {
        --tab-caption-width: 20%;
      }

      .install-instructions .sd-tab-set.tabs-os::before {
        content: "Operating System";
      }

      .install-instructions .sd-tab-set.tabs-package-manager::before {
        content: "Package Manager";
      }
    }
  </style>

.. div:: install-instructions

  .. tab-set::
    :class: tabs-os

    .. tab-item:: Windows
      :class-label: tab-4

      .. tab-set::
        :class: tabs-package-manager

        .. tab-item:: pip
          :class-label: tab-6
          :sync: package-manager-pip

          Install the 64-bit version of Python (|PythonMinVersion| or later), for instance from the
          `official website <https://www.python.org/downloads/windows/>`__.

          Now create a virtual environment (venv_) and install the required python packages:

          .. prompt:: powershell

            python -m venv sklearn-dev
            sklearn-dev\Scripts\activate  # activate
            pip install wheel numpy scipy cython meson-python ninja

          Also install the development dependencies:

          .. prompt:: powershell

            pip install pytest pytest-cov ruff==0.11.2 mypy numpydoc

          Additionally, you need to install a compiler with OpenMP_ support.
          Download the `Build Tools for Visual Studio installer <https://aka.ms/vs/17/release/vs_buildtools.exe>`_
          and run the downloaded `vs_buildtools.exe` file. During the installation you will
          need to make sure you select "Desktop development with C++", similarly to this
          screenshot:

          .. image::
            ../images/visual-studio-build-tools-selection.png

        .. tab-item:: conda
          :class-label: tab-6
          :sync: package-manager-conda

          Install a recent version of Python (|PythonMinVersion| or later) for instance
          using conda-forge_. Conda-forge provides a conda-based distribution of
          Python and the most popular scientific libraries.
          Create a new conda environment with the required python packages:

          .. prompt:: powershell

            conda create -n sklearn-dev -c conda-forge python numpy scipy cython meson-python ninja

          It is not always necessary but it is safer to open a new prompt before
          activating the newly created conda environment:

          .. prompt:: powershell

            conda activate sklearn-dev

          Also install the development dependencies in your environment:

          .. prompt:: powershell

            conda install -c conda-forge pytest pytest-cov ruff==0.11.2 mypy numpydoc

          Additionally, you need to install a compiler with OpenMP_ support.
          Download the `Build Tools for Visual Studio installer <https://aka.ms/vs/17/release/vs_buildtools.exe>`_
          and run the downloaded `vs_buildtools.exe` file. During the installation you will
          need to make sure you select "Desktop development with C++", similarly to this
          screenshot:

          .. image::
            ../images/visual-studio-build-tools-selection.png

    .. tab-item:: MacOS
      :class-label: tab-4

      .. tab-set::
        :class: tabs-package-manager

        .. tab-item:: pip
          :class-label: tab-6
          :sync: package-manager-pip

          Install Python 3 (|PythonMinVersion| or later) using Homebrew_
          (`brew install python`) or by manually installing the package from the
          `official website <https://www.python.org/downloads/macos/>`__.

          Now create a virtual environment (venv_) and install the required python packages:

          .. prompt:: bash

            python -m venv sklearn-dev
            source sklearn-dev/bin/activate  # activate
            pip install wheel numpy scipy cython meson-python ninja

          Also install the development dependencies:

          .. prompt:: bash

            pip install pytest pytest-cov ruff==0.11.2 mypy numpydoc

          The default C compiler on macOS, Apple clang (confusingly aliased as
          `/usr/bin/gcc`), does not directly support OpenMP, so you additionally
          need to enable OpenMP support.

          First install the macOS command line tools:

          .. prompt:: bash $

              xcode-select --install

          Install the LLVM OpenMP library with Homebrew_:

          .. prompt:: bash $

              brew install libomp

          Remove any existing scikit-learn installations and meson builds to avoid conflicts.
          You can use the provided `Makefile` for this by simply calling:

          .. prompt:: bash $

              make clean

        .. tab-item:: conda
          :class-label: tab-6
          :sync: package-manager-conda

          Install a recent version of Python (|PythonMinVersion| or later) for instance
          using conda-forge_. Conda-forge provides a conda-based distribution of
          Python and the most popular scientific libraries.
          Create a new conda environment with the required python packages:

          .. prompt:: bash $

            conda create -n sklearn-dev -c conda-forge python numpy scipy cython meson-python ninja

          It is not always necessary but it is safer to open a new prompt before
          activating the newly created conda environment:

          .. prompt:: bash $

            conda activate sklearn-dev

          Also install the development dependencies in your environment:

          .. prompt:: bash $

            conda install -c conda-forge pytest pytest-cov ruff==0.11.2 mypy numpydoc

          The default C compiler on macOS, Apple clang (confusingly aliased as
          `/usr/bin/gcc`), does not directly support OpenMP, so you need to
          install the ``compilers`` meta-package from the conda-forge channel,
          which provides OpenMP-enabled C/C++ compilers based on the llvm toolchain.

          First install the macOS command line tools:

          .. prompt:: bash $

              xcode-select --install

          Make sure you activated the `sklearn-dev` environment and install the following packages:

          .. prompt:: bash $

              conda install -c conda-forge joblib threadpoolctl pytest compilers llvm-openmp

          Remove any existing scikit-learn installations and meson builds to avoid conflicts.
          You can use the provided `Makefile` for this by simply calling:

          .. prompt:: bash $

              make clean

          .. note::

            If you get any conflicting dependency error message, try commenting out
            any custom conda configuration in the ``$HOME/.condarc`` file. In
            particular the ``channel_priority: strict`` directive is known to cause
            problems for this setup.

            You can check that the custom compilers are properly installed from conda
            forge using the following command:

            .. prompt:: bash $

                conda list

            which should include ``compilers`` and ``llvm-openmp``.

            The compilers meta-package will automatically set custom environment
            variables:

            .. prompt:: bash $

                echo $CC
                echo $CXX
                echo $CFLAGS
                echo $CXXFLAGS
                echo $LDFLAGS

            They point to files and folders from your ``sklearn-dev`` conda environment
            (in particular in the bin/, include/ and lib/ subfolders). For instance
            ``-L/path/to/conda/envs/sklearn-dev/lib`` should appear in ``LDFLAGS``.

            When installing scikit-learn in the next step, you should see the
            compiled extension being built with the clang and clang++ compilers installed by
            conda with the ``-fopenmp`` command line flag in the log.

    .. tab-item:: Linux
      :class-label: tab-4

      .. tab-set::
        :class: tabs-package-manager

        .. tab-item:: pip
          :class-label: tab-6
          :sync: package-manager-pip

          scikit-learn requires a recent version of Python (|PythonMinVersion| or later).
          Python 3 is usually installed by default on most Linux distributions. To
          check if you have it installed, try:

          .. prompt:: bash

            python3 --version

          If you don't have Python 3 installed, please install `python3` from your
          distribution's package manager.

          Now create a virtual environment (venv_) and install the required python packages:

          .. prompt:: bash

            python -m venv sklearn-dev
            source sklearn-dev/bin/activate  # activate
            pip install wheel numpy scipy cython meson-python ninja

          .. dropdown:: If no isolated environment is used (not recommended!)

            If, for some reason, you are not using an isolated environment (which is not recommended),
            you need to use `pip3` instead of `pip`. You can check if it is installed with:

            .. prompt:: bash

              pip3 --version

            and if not, install `python3-pip` from your distribution's package manager.

            Next, install the required python packages with:

            .. prompt:: bash

              pip3 install wheel numpy scipy cython meson-python ninja

            cython and the pre-compiled wheels for the runtime dependencies (numpy, scipy
            and joblib) should then automatically be installed in
            ``$HOME/.local/lib/pythonX.Y/site-packages``. In this case,
            ``pip3`` also needs to be used instead of ``pip`` in the following steps and
            when `installing scikit-learn <pip_build_>`_.

          Also install the development dependencies:

          .. prompt:: bash

            pip install pytest pytest-cov ruff==0.11.2 mypy numpydoc

          Additionally, you need to have a working C/C++ compiler with OpenMP support.

          Install the build dependencies for **Debian-based operating systems, e.g. Ubuntu**:

          .. prompt:: bash $

              sudo apt-get install build-essential python3-dev

          When precompiled wheels of the runtime dependencies are not available for your
          architecture (e.g. **ARM**), you can install the system versions:

          .. prompt:: bash $

              sudo apt-get install cython3 python3-numpy python3-scipy

          On **Red Hat and clones (e.g. CentOS)**, install the dependencies using:

          .. prompt:: bash $

              sudo yum -y install gcc gcc-c++ python3-devel numpy scipy

        .. tab-item:: conda
          :class-label: tab-6
          :sync: package-manager-conda

          scikit-learn requires a recent version of Python (|PythonMinVersion| or later).
          Python 3 is usually installed by default on most Linux distributions. To
          check if you have it installed, try:

          .. prompt:: bash

            python3 --version

          If you don't have Python |PythonMinVersion| or later, please install a
          recent version for instance using conda-forge_. Conda-forge provides a
          conda-based distribution of Python and the most popular scientific libraries.
          Create a new conda environment with the required python packages:

          .. prompt:: bash $

            conda create -n sklearn-dev -c conda-forge python numpy scipy cython meson-python ninja

          It is not always necessary but it is safer to open a new prompt before
          activating the newly created conda environment:

          .. prompt:: bash $

            conda activate sklearn-dev

          Also install the development dependencies in your environment:

          .. prompt:: bash $

            conda install -c conda-forge pytest pytest-cov ruff==0.11.2 mypy numpydoc

          Additionally, you need to have the scikit-learn Python development headers
          and a working C/C++ compiler with OpenMP support:

          .. prompt:: bash $

              conda install -c conda-forge joblib threadpoolctl compilers

          Remove any existing scikit-learn installations and meson builds to avoid conflicts.
          You can use the provided `Makefile` for this by simply calling:

          .. prompt:: bash $

              make clean

Install editable version of scikit-learn
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Within your venv or conda `sklearn-dev` environment, build the project with pip:

.. prompt:: bash $

  pip install --editable . \
    --verbose --no-build-isolation \
    --config-settings editable-verbose=true

.. note::

  `--config-settings editable-verbose=true` is optional but recommended
  to avoid surprises when you import `sklearn`. `meson-python` implements
  editable installs by rebuilding `sklearn` when executing `import sklearn`.
  With the recommended setting you will see a message when this happens,
  rather than potentially waiting without feedback and wondering
  what is taking so long. Bonus: this means you only have to run the `pip
  install` command once, `sklearn` will automatically be rebuilt when
  importing `sklearn`.

  Note that `--config-settings` is only supported in `pip` version 23.1 or
  later. To upgrade `pip` to a compatible version, run `pip install -U pip`.

To check your installation, make sure that the installed scikit-learn has a
version number ending with `.dev0`:

  .. prompt:: bash $

    python -c "import sklearn; sklearn.show_versions()"

You should now have a working installation of scikit-learn and your git repository
properly configured.
Please refer to the :ref:`developers_guide` and :ref:`pytest_tips` to run
some tests to verify your installation.

.. _OpenMP: https://en.wikipedia.org/wiki/OpenMP
.. _Cython: https://cython.org
.. _meson-python: https://mesonbuild.com/meson-python
.. _Ninja: https://ninja-build.org/
.. _NumPy: https://numpy.org
.. _SciPy: https://www.scipy.org
.. _Homebrew: https://brew.sh
.. _venv: https://docs.python.org/3/tutorial/venv.html
.. _conda environment: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html
.. _conda-forge: https://conda-forge.org/download/
