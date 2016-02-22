
.. _advanced-installation:

===================================
Advanced installation instructions
===================================

There are different ways to get scikit-learn installed:

  * Install the version of scikit-learn provided by your
    :ref:`operating system or Python distribution <install_by_distribution>`.
    This is the quickest option for those who have operating systems that
    distribute scikit-learn.

  * :ref:`Install an official release <install_official_release>`. This
    is the best approach for users who want a stable version number
    and aren't concerned about running a slightly older version of
    scikit-learn.

  * :ref:`Install the latest development version
    <install_bleeding_edge>`. This is best for users who want the
    latest-and-greatest features and aren't afraid of running
    brand-new code.

.. note::

    If you wish to contribute to the project, you need to
    :ref:`install the latest development version<install_bleeding_edge>`.


.. _install_official_release:

Installing an official release
==============================

Scikit-learn requires:

- Python (>= 2.6 or >= 3.3),
- NumPy (>= 1.6.1),
- SciPy (>= 0.9).


Mac OSX
-------

Scikit-learn and its dependencies are all available as wheel packages for OSX::

    pip install -U numpy scipy scikit-learn


Linux
-----

At this time scikit-learn does not provide official binary packages for Linux
so you have to build from source if you want the lastest version.
If you don't need the newest version, consider using your package manager to
install scikit-learn. it is usually the easiest way, but might not provide the
newest version.

installing build dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

installing from source requires you to have installed the scikit-learn runtime
dependencies, python development headers and a working c/c++ compiler.
under debian-based operating systems, which include ubuntu, if you have
python 2 you can install all these requirements by issuing::

    sudo apt-get install build-essential python-dev python-setuptools \
                         python-numpy python-scipy \
                         libatlas-dev libatlas3gf-base

if you have python 3::

    sudo apt-get install build-essential python3-dev python3-setuptools \
                         python3-numpy python3-scipy \
                         libatlas-dev libatlas3gf-base

on recent debian and ubuntu (e.g. ubuntu 13.04 or later) make sure that atlas
is used to provide the implementation of the blas and lapack linear algebra
routines::

    sudo update-alternatives --set libblas.so.3 \
        /usr/lib/atlas-base/atlas/libblas.so.3
    sudo update-alternatives --set liblapack.so.3 \
        /usr/lib/atlas-base/atlas/liblapack.so.3

.. note::

    in order to build the documentation and run the example code contains in
    this documentation you will need matplotlib::

        sudo apt-get install python-matplotlib

.. note::

    the above installs the atlas implementation of blas
    (the basic linear algebra subprograms library).
    ubuntu 11.10 and later, and recent (testing) versions of debian,
    offer an alternative implementation called openblas.

    using openblas can give speedups in some scikit-learn modules,
    but can freeze joblib/multiprocessing prior to openblas version 0.2.8-4,
    so using it is not recommended unless you know what you're doing.

    if you do want to use openblas, then replacing atlas only requires a couple
    of commands. atlas has to be removed, otherwise numpy may not work::

        sudo apt-get remove libatlas3gf-base libatlas-dev
        sudo apt-get install libopenblas-dev

        sudo update-alternatives  --set libblas.so.3 \
            /usr/lib/openblas-base/libopenblas.so.0
        sudo update-alternatives --set liblapack.so.3 \
            /usr/lib/lapack/liblapack.so.3

on red hat and clones (e.g. centos), install the dependencies using::

    sudo yum -y install gcc gcc-c++ numpy python-devel scipy


building scikit-learn with pip
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

this is usually the fastest way to install or upgrade to the latest stable
release::

    pip install --user --install-option="--prefix=" -u scikit-learn

the ``--user`` flag asks pip to install scikit-learn in the ``$home/.local``
folder therefore not requiring root permission. this flag should make pip
ignore any old version of scikit-learn previously installed on the system while
benefiting from system packages for numpy and scipy. those dependencies can
be long and complex to build correctly from source.

the ``--install-option="--prefix="`` flag is only required if python has a
``distutils.cfg`` configuration with a predefined ``prefix=`` entry.


from source package
~~~~~~~~~~~~~~~~~~~

download the source package from 
`pypi <https://pypi.python.org/pypi/scikit-learn>`_,
, unpack the sources and cd into the source directory.

this packages uses distutils, which is the default way of installing
python modules. the install command is::

    python setup.py install

or alternatively (also from within the scikit-learn source folder)::

    pip install .

.. warning::

   packages installed with the ``python setup.py install`` command cannot
   be uninstalled nor upgraded by ``pip`` later. to properly uninstall
   scikit-learn in that case it is necessary to delete the ``sklearn`` folder
   from your python ``site-packages`` directory.


windows
-------

first, you need to install `numpy <http://www.numpy.org/>`_ and `scipy
<http://www.scipy.org/>`_ from their own official installers.

wheel packages (.whl files) for scikit-learn from `pypi
<https://pypi.python.org/pypi/scikit-learn/>`_ can be installed with the `pip
<https://pip.readthedocs.org/en/stable/installing/>`_ utility.
open a console and type the following to install or upgrade scikit-learn to the
latest stable release::

    pip install -u scikit-learn

if there are no binary packages matching your python, version you might
to try to install scikit-learn and its dependencies from `christoph gohlke
unofficial windows installers
<http://www.lfd.uci.edu/~gohlke/pythonlibs/#scikit-learn>`_
or from a :ref:`python distribution <install_by_distribution>` instead.


.. _install_by_distribution:

third party distributions of scikit-learn
=========================================

some third-party distributions are now providing versions of
scikit-learn integrated with their package-management systems.

these can make installation and upgrading much easier for users since
the integration includes the ability to automatically install
dependencies (numpy, scipy) that scikit-learn requires.

the following is an incomplete list of python and os distributions
that provide their own version of scikit-learn.


macports for mac osx
--------------------

the macports package is named ``py<xy>-scikits-learn``,
where ``xy`` denotes the python version.
it can be installed by typing the following
command::

    sudo port install py26-scikit-learn

or::

    sudo port install py27-scikit-learn


arch linux
----------

arch linux's package is provided through the `official repositories
<https://www.archlinux.org/packages/?q=scikit-learn>`_ as
``python-scikit-learn`` for python 3 and ``python2-scikit-learn`` for python 2.
it can be installed by typing the following command:

.. code-block:: none

     # pacman -s python-scikit-learn

or:

.. code-block:: none

     # pacman -s python2-scikit-learn

depending on the version of python you use.


netbsd
------

scikit-learn is available via `pkgsrc-wip <http://pkgsrc-wip.sourceforge.net/>`_:

    http://pkgsrc.se/wip/py-scikit_learn

fedora
------

the fedora package is called ``python-scikit-learn`` for the python 2 version
and ``python3-scikit-learn`` for the python 3 version. both versions can
be installed using ``yum``::

    $ sudo yum install python-scikit-learn

or::

    $ sudo yum install python3-scikit-learn


building on windows
-------------------

to build scikit-learn on windows you need a working c/c++ compiler in
addition to numpy, scipy and setuptools.

picking the right compiler depends on the version of python (2 or 3)
and the architecture of the python interpreter, 32-bit or 64-bit.
you can check the python version by running the following in ``cmd`` or
``powershell`` console::

    python --version

and the architecture with::

    python -c "import struct; print(struct.calcsize('p') * 8)"

the above commands assume that you have the python installation folder in your
path environment variable.


32-bit python
-------------

for 32-bit python it is possible use the standalone installers for
`microsoft visual c++ express 2008 <http://go.microsoft.com/?linkid=7729279>`_
for python 2 or microsoft visual c++ express 2010 for python 3.

once installed you should be able to build scikit-learn without any
particular configuration by running the following command in the scikit-learn
folder::

   python setup.py install


64-bit python
-------------

for the 64-bit architecture, you either need the full visual studio or
the free windows sdks that can be downloaded from the links below.

the windows sdks include the msvc compilers both for 32 and 64-bit
architectures. they come as a ``grmsdkx_en_dvd.iso`` file that can be mounted
as a new drive with a ``setup.exe`` installer in it.

- for python 2 you need sdk **v7.0**: `ms windows sdk for windows 7 and .net
  framework 3.5 sp1
  <http://www.microsoft.com/en-us/download/details.aspx?id=18950>`_

- for python 3 you need sdk **v7.1**: `ms windows sdk for windows 7 and .net
  framework 4
  <https://www.microsoft.com/en-us/download/details.aspx?id=8442>`_

both sdks can be installed in parallel on the same host. to use the windows
sdks, you need to setup the environment of a ``cmd`` console launched with the
following flags (at least for sdk v7.0)::

    cmd /e:on /v:on /k

then configure the build environment with::

    set distutils_use_sdk=1
    set mssdk=1
    "c:\program files\microsoft sdks\windows\v7.0\setup\windowssdkver.exe" -q -version:v7.0
    "c:\program files\microsoft sdks\windows\v7.0\bin\setenv.cmd" /x64 /release

finally you can build scikit-learn in the same ``cmd`` console::

    python setup.py install

replace ``v7.0`` by the ``v7.1`` in the above commands to do the same for
python 3 instead of python 2.

replace ``/x64`` by ``/x86``  to build for 32-bit python instead of 64-bit
python.


building binary packages and installers
---------------------------------------

the ``.whl`` package and ``.exe`` installers can be built with::

    pip install wheel
    python setup.py bdist_wheel bdist_wininst -b doc/logos/scikit-learn-logo.bmp

the resulting packages are generated in the ``dist/`` folder.


using an alternative compiler
-----------------------------

it is possible to use `mingw <http://www.mingw.org>`_ (a port of gcc to windows
os) as an alternative to msvc for 32-bit python. not that extensions built with
mingw32 can be redistributed as reusable packages as they depend on gcc runtime
libraries typically not installed on end-users environment.

to force the use of a particular compiler, pass the ``--compiler`` flag to the
build step::

    python setup.py build --compiler=my_compiler install

where ``my_compiler`` should be one of ``mingw32`` or ``msvc``.


.. _install_bleeding_edge:

bleeding edge
=============

see section :ref:`git_repo` on how to get the development version. then follow
the previous instructions to build from source depending on your platform. 
You will also require Cython >=0.23 in order to build the development version.


.. _testing:

testing
=======

testing scikit-learn once installed
-----------------------------------

testing requires having the `nose
<https://nose.readthedocs.org/en/latest/>`_ library. after
installation, the package can be tested by executing *from outside* the
source directory::

    $ nosetests -v sklearn

under windows, it is recommended to use the following command (adjust the path
to the ``python.exe`` program) as using the ``nosetests.exe`` program can badly
interact with tests that use ``multiprocessing``::

    c:\python34\python.exe -c "import nose; nose.main()" -v sklearn

this should give you a lot of output (and some warnings) but
eventually should finish with a message similar to::

    ran 3246 tests in 260.618s
    ok (skip=20)

otherwise, please consider posting an issue into the `bug tracker
<https://github.com/scikit-learn/scikit-learn/issues>`_ or to the
:ref:`mailing_lists` including the traceback of the individual failures
and errors. please include your operation system, your version of numpy, scipy
and scikit-learn, and how you installed scikit-learn.


testing scikit-learn from within the source folder
--------------------------------------------------

scikit-learn can also be tested without having the package
installed. for this you must compile the sources inplace from the
source directory::

    python setup.py build_ext --inplace

test can now be run using nosetests::

    nosetests -v sklearn/

this is automated by the commands::

    make in

and::

    make test


you can also install a symlink named ``site-packages/scikit-learn.egg-link``
to the development folder of scikit-learn with::

    pip install --editable .
