.. _installation-instructions:

=======================
Installing scikit-learn
=======================

There are different ways to install scikit-learn:

  * :ref:`Install the latest official release <install_official_release>`. This
    is the best approach for most users. It will provide a stable version
    and pre-built packages are available for most platforms.

  * Install the version of scikit-learn provided by your
    :ref:`operating system or Python distribution <install_by_distribution>`.
    This is a quick option for those who have operating systems or Python
    distributions that distribute scikit-learn.
    It might not provide the latest release version.

  * :ref:`Building the package from source
    <install_bleeding_edge>`. This is best for users who want the
    latest-and-greatest features and aren't afraid of running
    brand-new code. This is also needed for users who wish to contribute to the
    project.


.. _install_official_release:

Installing the latest release
=============================

.. This quickstart installation is a hack of the awesome
   https://spacy.io/usage/#quickstart page.
   See the original javascript implementation
   https://github.com/ines/quickstart


.. raw:: html

  <div class="install">
       <strong>Operating System</strong>
          <input type="radio" name="os" id="quickstart-win" checked>
          <label for="quickstart-win">Windows</label>
          <input type="radio" name="os" id="quickstart-mac">
          <label for="quickstart-mac">macOS</label>
          <input type="radio" name="os" id="quickstart-lin">
          <label for="quickstart-lin">Linux</label><br />
       <strong>Packager</strong>
          <input type="radio" name="packager" id="quickstart-pip" checked>
          <label for="quickstart-pip">pip</label>
          <input type="radio" name="packager" id="quickstart-conda">
          <label for="quickstart-conda">conda</label><br />
          <input type="checkbox" name="config" id="quickstart-venv">
          <label for="quickstart-venv"></label>
       </span>

.. raw:: html

       <div>
         <span class="sk-expandable" data-packager="pip" data-os="windows">Install the 64bit version of Python 3, for instance from <a href="https://www.python.org/">https://www.python.org</a>.</span
         ><span class="sk-expandable" data-packager="pip" data-os="mac">Install Python 3 using <a href="https://brew.sh/">homebrew</a> (<code>brew install python</code>) or by manually installing the package from <a href="https://www.python.org">https://www.python.org</a>.</span
         ><span class="sk-expandable" data-packager="pip" data-os="linux">Install python3 and python3-pip using the package manager of the Linux Distribution.</span
         ><span class="sk-expandable" data-packager="conda"
            >Install conda using the <a href="https://docs.conda.io/projects/conda/en/latest/user-guide/install/">Anaconda or miniconda</a>
             installers or the <a href="https://https://github.com/conda-forge/miniforge#miniforge">miniforge</a> installers
             (no administrator permission required for any of those).</span>
       </div>

Then run:

.. raw:: html

       <div class="highlight"><pre><code
        ><span class="sk-expandable" data-packager="pip" data-os="linux" data-venv="">python3 -m venv sklearn-venv</span
        ><span class="sk-expandable" data-packager="pip" data-os="windows" data-venv="">python -m venv sklearn-venv</span
        ><span class="sk-expandable" data-packager="pip" data-os="mac" data-venv="">python -m venv sklearn-venv</span
        ><span class="sk-expandable" data-packager="pip" data-os="linux" data-venv="">source sklearn-venv/bin/activate</span
        ><span class="sk-expandable" data-packager="pip" data-os="mac" data-venv="">source sklearn-venv/bin/activate</span
        ><span class="sk-expandable" data-packager="pip" data-os="windows" data-venv="">sklearn-venv\Scripts\activate</span
        ><span class="sk-expandable" data-packager="pip" data-venv="">pip install -U scikit-learn</span
        ><span class="sk-expandable" data-packager="pip" data-os="mac" data-venv="no">pip install -U scikit-learn</span
        ><span class="sk-expandable" data-packager="pip" data-os="windows" data-venv="no">pip install -U scikit-learn</span
        ><span class="sk-expandable" data-packager="pip" data-os="linux" data-venv="no">pip3 install -U scikit-learn</span
        ><span class="sk-expandable" data-packager="conda" data-venv="">conda create -n sklearn-env</span
        ><span class="sk-expandable" data-packager="conda" data-venv="">conda activate sklearn-env</span
        ><span class="sk-expandable" data-packager="conda">conda install -c conda-forge scikit-learn </span
       ></code></pre></div>

In order to check your installation you can use

.. raw:: html

   <div class="highlight"><pre><code
      ><span class="sk-expandable" data-packager="pip" data-os="linux" data-venv="no">python3 -m pip show scikit-learn  # to see which version and where scikit-learn is installed</span
      ><span class="sk-expandable" data-packager="pip" data-os="linux" data-venv="no">python3 -m pip freeze  # to see all packages installed in the active virtualenv</span
      ><span class="sk-expandable" data-packager="pip" data-os="linux" data-venv="no">python3 -c "import sklearn; sklearn.show_versions()"</span
      ><span class="sk-expandable" data-packager="pip" data-venv="">python -m pip show scikit-learn  # to see which version and where scikit-learn is installed</span
      ><span class="sk-expandable" data-packager="pip" data-venv="">python -m pip freeze  # to see all packages installed in the active virtualenv</span
      ><span class="sk-expandable" data-packager="pip" data-venv="">python -c "import sklearn; sklearn.show_versions()"</span
      ><span class="sk-expandable" data-packager="pip" data-os="windows" data-venv="no">python -m pip show scikit-learn  # to see which version and where scikit-learn is installed</span
      ><span class="sk-expandable" data-packager="pip" data-os="windows" data-venv="no">python -m pip freeze  # to see all packages installed in the active virtualenv</span
      ><span class="sk-expandable" data-packager="pip" data-os="windows" data-venv="no">python -c "import sklearn; sklearn.show_versions()"</span
      ><span class="sk-expandable" data-packager="pip" data-os="mac" data-venv="no">python -m pip show scikit-learn  # to see which version and where scikit-learn is installed</span
      ><span class="sk-expandable" data-packager="pip" data-os="mac" data-venv="no">python -m pip freeze  # to see all packages installed in the active virtualenv</span
      ><span class="sk-expandable" data-packager="pip" data-os="mac" data-venv="no">python -c "import sklearn; sklearn.show_versions()"</span
      ><span class="sk-expandable" data-packager="conda">conda list scikit-learn  # to see which scikit-learn version is installed</span
      ><span class="sk-expandable" data-packager="conda">conda list  # to see all packages installed in the active conda environment</span
      ><span class="sk-expandable" data-packager="conda">python -c "import sklearn; sklearn.show_versions()"</span
      ></code></pre></div>
  </div>

Note that in order to avoid potential conflicts with other packages it is
strongly recommended to use a `virtual environment (venv)
<https://docs.python.org/3/tutorial/venv.html>`_ or a `conda environment
<https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_.

Using such an isolated environment makes it possible to install a specific
version of scikit-learn with pip or conda and its dependencies independently of
any previously installed Python packages. In particular under Linux is it
discouraged to install pip packages alongside the packages managed by the
package manager of the distribution (apt, dnf, pacman...).

Note that you should always remember to activate the environment of your choice
prior to running any Python command whenever you start a new terminal session.

If you have not installed NumPy or SciPy yet, you can also install these using
conda or pip. When using pip, please ensure that *binary wheels* are used,
and NumPy and SciPy are not recompiled from source, which can happen when using
particular configurations of operating system and hardware (such as Linux on
a Raspberry Pi).


Scikit-learn plotting capabilities (i.e., functions start with "plot\_"
and classes end with "Display") require Matplotlib. The examples require
Matplotlib and some examples require scikit-image, pandas, or seaborn. The
minimum version of Scikit-learn dependencies are listed below along with its
purpose.

.. include:: min_dependency_table.rst

.. warning::

    Scikit-learn 0.20 was the last version to support Python 2.7 and Python 3.4.
    Scikit-learn 0.21 supported Python 3.5-3.7.
    Scikit-learn 0.22 supported Python 3.5-3.8.
    Scikit-learn now requires Python 3.6 or newer.


.. note::

   For installing on PyPy, PyPy3-v5.10+, Numpy 1.14.0+, and scipy 1.1.0+
   are required.

.. _install_on_apple_silicon_m1:

Installing on Apple Silicon M1 hardware
=======================================

The recently introduced `macos/arm64` platform (sometimes also known as
`macos/aarch64`) requires the open source community to upgrade the build
configuation and automation to properly support it.

At the time of writing (January 2021), the only way to get a working
installation of scikit-learn on this hardware is to install scikit-learn and its
dependencies from the conda-forge distribution, for instance using the miniforge
installers:

https://github.com/conda-forge/miniforge

The following issue tracks progress on making it possible to install
scikit-learn from PyPI with pip:

https://github.com/scikit-learn/scikit-learn/issues/19137


.. _install_by_distribution:

Third party distributions of scikit-learn
=========================================

Some third-party distributions provide versions of
scikit-learn integrated with their package-management systems.

These can make installation and upgrading much easier for users since
the integration includes the ability to automatically install
dependencies (numpy, scipy) that scikit-learn requires.

The following is an incomplete list of OS and python distributions
that provide their own version of scikit-learn.

Arch Linux
----------

Arch Linux's package is provided through the `official repositories
<https://www.archlinux.org/packages/?q=scikit-learn>`_ as
``python-scikit-learn`` for Python.
It can be installed by typing the following command:

.. prompt:: bash $

  sudo pacman -S python-scikit-learn


Debian/Ubuntu
-------------

The Debian/Ubuntu package is splitted in three different packages called
``python3-sklearn`` (python modules), ``python3-sklearn-lib`` (low-level
implementations and bindings), ``python3-sklearn-doc`` (documentation).
Only the Python 3 version is available in the Debian Buster (the more recent
Debian distribution).
Packages can be installed using ``apt-get``:

.. prompt:: bash $

  sudo apt-get install python3-sklearn python3-sklearn-lib python3-sklearn-doc


Fedora
------

The Fedora package is called ``python3-scikit-learn`` for the python 3 version,
the only one available in Fedora30.
It can be installed using ``dnf``:

.. prompt:: bash $

  sudo dnf install python3-scikit-learn


NetBSD
------

scikit-learn is available via `pkgsrc-wip
<http://pkgsrc-wip.sourceforge.net/>`_:

    http://pkgsrc.se/math/py-scikit-learn


MacPorts for Mac OSX
--------------------

The MacPorts package is named ``py<XY>-scikits-learn``,
where ``XY`` denotes the Python version.
It can be installed by typing the following
command:

.. prompt:: bash $

  sudo port install py36-scikit-learn


Anaconda and Enthought Deployment Manager for all supported platforms
---------------------------------------------------------------------

`Anaconda <https://www.anaconda.com/download>`_ and
`Enthought Deployment Manager <https://assets.enthought.com/downloads/>`_
both ship with scikit-learn in addition to a large set of scientific
python library for Windows, Mac OSX and Linux.

Anaconda offers scikit-learn as part of its free distribution.


Intel conda channel
-------------------

Intel maintains a dedicated conda channel that ships scikit-learn:

.. prompt:: bash $

  conda install -c intel scikit-learn

This version of scikit-learn comes with alternative solvers for some common
estimators. Those solvers come from the DAAL C++ library and are optimized for
multi-core Intel CPUs.

Note that those solvers are not enabled by default, please refer to the
`daal4py <https://intelpython.github.io/daal4py/sklearn.html>`_ documentation
for more details.

Compatibility with the standard scikit-learn solvers is checked by running the
full scikit-learn test suite via automated continuous integration as reported
on https://github.com/IntelPython/daal4py.


WinPython for Windows
-----------------------

The `WinPython <https://winpython.github.io/>`_ project distributes
scikit-learn as an additional plugin.


Troubleshooting
===============

.. _windows_longpath:

Error caused by file path length limit on Windows
-------------------------------------------------

It can happen that pip fails to install packages when reaching the default path
size limit of Windows if Python is installed in a nested location such as the
`AppData` folder structure under the user home directory, for instance::

    C:\Users\username>C:\Users\username\AppData\Local\Microsoft\WindowsApps\python.exe -m pip install scikit-learn
    Collecting scikit-learn
    ...
    Installing collected packages: scikit-learn
    ERROR: Could not install packages due to an EnvironmentError: [Errno 2] No such file or directory: 'C:\\Users\\username\\AppData\\Local\\Packages\\PythonSoftwareFoundation.Python.3.7_qbz5n2kfra8p0\\LocalCache\\local-packages\\Python37\\site-packages\\sklearn\\datasets\\tests\\data\\openml\\292\\api-v1-json-data-list-data_name-australian-limit-2-data_version-1-status-deactivated.json.gz'

In this case it is possible to lift that limit in the Windows registry by
using the ``regedit`` tool:

#. Type "regedit" in the Windows start menu to launch ``regedit``.

#. Go to the
   ``Computer\HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\FileSystem``
   key.

#. Edit the value of the ``LongPathsEnabled`` property of that key and set
   it to 1.

#. Reinstall scikit-learn (ignoring the previous broken installation):

.. prompt:: python $

    pip install --exists-action=i scikit-learn
