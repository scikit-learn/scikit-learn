.. _install_scikit_learn:

Installing scikit-learn
=======================

Scikit-learn can be installed from PyPI (using pip) or from conda-forge (using
conda). Both methods are fully supported and widely used. If you already use
conda, conda-forge provides optimized packages for scientific computing. If you
prefer pip, PyPI wheels are equally valid.

We recommend installing scikit-learn in an isolated environment to avoid
conflicts with other packages. Below are basic instructions for common operating
systems.

Installing the latest release
=============================

To install scikit-learn, follow the instructions for your operating system:

.. raw:: html

  <style>
    /* Show caption on large screens */
    @media screen and (min-width: 960px) {
      .install-instructions .sd-tab-set {
        --tab-caption-width: 20%;
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

          1. Download and install the 64-bit Python 3 from the
             `Windows official website <https://www.python.org/downloads/windows/>`_.

          2. Create a virtual environment (recommended):

             .. prompt:: powershell

                python -m venv sklearn-env
                sklearn-env\Scripts\activate

          3. Install scikit-learn using pip:

             .. prompt:: powershell

                pip install -U scikit-learn

          4. Verify the installation:

             .. prompt:: powershell

                python -m pip show scikit-learn
                python -c "import sklearn; sklearn.show_versions()"

        .. tab-item:: conda
          :class-label: tab-6

          1. Install Conda using the
             `Miniforge installers <https://github.com/conda-forge/miniforge#download>`_.

          2. Create a new Conda environment and install scikit-learn:

             .. prompt:: powershell

                conda create -n sklearn-env -c conda-forge scikit-learn
                conda activate sklearn-env

          3. Verify the installation:

             .. prompt:: powershell

                conda list scikit-learn
                python -c "import sklearn; sklearn.show_versions()"

    .. tab-item:: macOS
      :class-label: tab-4

      .. tab-set::
        :class: tabs-package-manager

        .. tab-item:: pip
          :class-label: tab-6

          1. Install Python 3 using `Homebrew <https://brew.sh>`_ or from the
             `macOS official website <https://www.python.org/downloads/macos/>`_.

          2. Create a virtual environment (recommended):

             .. prompt:: bash

                python3 -m venv sklearn-env
                source sklearn-env/bin/activate

          3. Install scikit-learn using pip:

             .. prompt:: bash

                pip install -U scikit-learn

          4. Verify the installation:

             .. prompt:: bash

                python3 -m pip show scikit-learn
                python3 -c "import sklearn; sklearn.show_versions()"

        .. tab-item:: conda
          :class-label: tab-6

          1. Install Conda using the
             `Miniforge installers <https://github.com/conda-forge/miniforge#download>`_.

          2. Create a new Conda environment and install scikit-learn:

             .. prompt:: bash

                conda create -n sklearn-env -c conda-forge scikit-learn
                conda activate sklearn-env

          3. Verify the installation:

             .. prompt:: bash

                conda list scikit-learn
                python3 -c "import sklearn; sklearn.show_versions()"

    .. tab-item:: Linux
      :class-label: tab-4

      .. tab-set::
        :class: tabs-package-manager

        .. tab-item:: pip
          :class-label: tab-6

          1. Ensure Python 3 and pip are installed:

             .. prompt:: bash

                python3 --version
                pip3 --version

          2. Create a virtual environment (recommended):

             .. prompt:: bash

                python3 -m venv sklearn-env
                source sklearn-env/bin/activate

          3. Install scikit-learn using pip:

             .. prompt:: bash

                pip3 install -U scikit-learn

          4. Verify the installation:

             .. prompt:: bash

                python3 -m pip show scikit-learn
                python3 -c "import sklearn; sklearn.show_versions()"

        .. tab-item:: conda
          :class-label: tab-6

          1. Install Conda using the
             `Miniforge installers <https://github.com/conda-forge/miniforge#download>`_.

          2. Create a new Conda environment and install scikit-learn:

             .. prompt:: bash

                conda create -n sklearn-env -c conda-forge scikit-learn
                conda activate sklearn-env

          3. Verify the installation:

             .. prompt:: bash

                conda list scikit-learn
                python3 -c "import sklearn; sklearn.show_versions()"

---

**Important**: Using an isolated environment such as ``pip venv`` or ``conda``
ensures that scikit-learn and its dependencies are installed independently of
other Python packages. Always activate the environment before running Python
commands to avoid conflicts.

For more details on Python packaging and virtual environments, refer to the
`PyPackaging Native documentation <https://py-pkgs.org>`_.


Third-party distributions of scikit-learn
=========================================

Several third-party distributions include scikit-learn in their package
managers. This can simplify installation and upgrading scikit-learn, as
dependencies (NumPy, SciPy, etc.) are often handled automatically. Below is an
incomplete list:

- **Alpine Linux** (`py3-scikit-learn docs
  <https://pkgs.alpinelinux.org/packages?name=py3-scikit-learn>`_):
  ``sudo apk add py3-scikit-learn``
- **Arch Linux** (`python-scikit-learn
  <https://archlinux.org/packages/?q=scikit-learn>`_):
  ``sudo pacman -S python-scikit-learn``
- **Debian/Ubuntu** (`python3-sklearn
  <https://packages.debian.org/search?keywords=python3-sklearn>`_):
  ``sudo apt-get install python3-sklearn python3-sklearn-lib python-sklearn-doc``
- **Fedora** (`python3-scikit-learn
  <https://apps.fedoraproject.org/packages/python3-scikit-learn>`_):
  ``sudo dnf install python3-scikit-learn``
- **NetBSD**: via `pkgsrc-wip <https://pkgsrc.se/math/py-scikit-learn>`_
- **macOS (MacPorts)** (`py-scikits-learn
  <https://ports.macports.org/port/py-scikits-learn/>`_):
  ``sudo port install py39-scikit-learn``
- **Intel Extension for scikit-learn** (`docs
  <https://intel.github.io/scikit-learn-intelex>`_): install via
  ``pip install scikit-learn-intelex`` or
  ``conda install scikit-learn-intelex``
- **WinPython for Windows** (`homepage
  <https://winpython.github.io/>`_)

.. note::
   Third-party repositories may not always provide the latest release. Consult
   each distributor’s documentation for more details.



Troubleshooting
===============

If you encounter unexpected failures when installing scikit-learn, you may submit
an issue to the `issue tracker <https://github.com/scikit-learn/scikit-learn/issues>`_.
Before that, please check these common issues:

Error caused by file path length limit on Windows
-------------------------------------------------

Pip can fail to install packages when hitting the default path-size limit on
Windows, especially if Python is installed under ``AppData``. For example::

   C:\Users\username>...
   ERROR: Could not install packages due to an OSError: [Errno 2] No such file
   or directory:
   'C:\\Users\\username\\AppData\\Local\\...\\sklearn\\datasets\\tests\\data\\...'

To fix this, enable long paths via the Windows registry:

1. Type ``regedit`` in the Start menu.
2. Go to ``Computer\HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\FileSystem``.
3. Edit the ``LongPathsEnabled`` property to ``1``.

Then reinstall scikit-learn (ignoring the prior broken install):

.. prompt:: powershell

   pip install --exists-action=i scikit-learn
