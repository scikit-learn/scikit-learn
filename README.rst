.. -*- mode: rst -*-

.. raw:: html

   <h1 align="center" style="margin-top: 0.5rem; margin-bottom: 0.5rem; font-size: 4rem;">
     Simple and efficient tools for predictive data analysis in Python
   </h1>

|GitHubActions| |CircleCI| |Codecov| |Nightly wheels| |Ruff| |PyPI| |PythonVersion| |DOI| |Benchmark|

.. |GitHubActions| image:: https://github.com/scikit-learn/scikit-learn/actions/workflows/unit-tests.yml/badge.svg?
   :target: https://github.com/scikit-learn/scikit-learn/actions/workflows/unit-tests.yml?query=branch%3Amain
   :alt: GitHub Actions: unit tests

.. |CircleCI| image:: https://circleci.com/gh/scikit-learn/scikit-learn/tree/main.svg?style=shield
   :target: https://circleci.com/gh/scikit-learn/scikit-learn
   :alt: CircleCI

.. |Codecov| image:: https://codecov.io/gh/scikit-learn/scikit-learn/branch/main/graph/badge.svg?token=Pk8G9gg3y9
   :target: https://codecov.io/gh/scikit-learn/scikit-learn
   :alt: Codecov

.. |Nightly wheels| image:: https://github.com/scikit-learn/scikit-learn/actions/workflows/wheels.yml/badge.svg?event=schedule
   :target: https://github.com/scikit-learn/scikit-learn/actions?query=workflow%3A%22Wheel+builder%22+event%3Aschedule
   :alt: Nightly wheels

.. |Ruff| image:: https://img.shields.io/badge/code%20style-ruff-000000.svg?
   :target: https://github.com/astral-sh/ruff
   :alt: Code style: Ruff

.. |PyPI| image:: https://img.shields.io/pypi/v/scikit-learn?color=ff8c00&label=PyPI&logo=pypi
   :target: https://pypi.org/project/scikit-learn
   :alt: PyPI version

.. |PythonVersion| image:: https://img.shields.io/pypi/pyversions/scikit-learn.svg?color=ff8c00
   :target: https://pypi.org/project/scikit-learn/
   :alt: Supported Python versions

.. |DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.17880109.svg?
   :target: https://zenodo.org/badge/latestdoi/21369/scikit-learn/scikit-learn
   :alt: DOI

.. |Benchmark| image:: https://img.shields.io/badge/Benchmarked%20by-asv-0a7cbc
   :target: https://scikit-learn.org/scikit-learn-benchmarks
   :alt: Benchmarked by asv

.. |PythonMinVersion| replace:: 3.11
.. |NumPyMinVersion| replace:: 1.24.1
.. |SciPyMinVersion| replace:: 1.10.0
.. |JoblibMinVersion| replace:: 1.4.0
.. |NarwhalsMinVersion| replace:: 2.0.1
.. |ThreadpoolctlMinVersion| replace:: 3.5.0
.. |MatplotlibMinVersion| replace:: 3.6.1
.. |Scikit-ImageMinVersion| replace:: 0.22.0
.. |PandasMinVersion| replace:: 1.5.0
.. |SeabornMinVersion| replace:: 0.13.0
.. |PytestMinVersion| replace:: 7.1.2
.. |PlotlyMinVersion| replace:: 5.22.0

**scikit-learn** is a Python module for machine learning built on top of
SciPy and is distributed under the 3-Clause BSD license.

The project was started in 2007 by David Cournapeau as a Google Summer
of Code project, and since then many volunteers have contributed. See
the `About us <https://scikit-learn.org/dev/about.html#authors>`__ page
for a list of core contributors.

It is currently maintained by a team of volunteers.

Website: https://scikit-learn.org

Installation
------------

Dependencies
~~~~~~~~~~~~

scikit-learn requires:

- Python (>= |PythonMinVersion|)
- NumPy (>= |NumPyMinVersion|)
- SciPy (>= |SciPyMinVersion|)
- Narwhals (>= |NarwhalsMinVersion|)
- joblib (>= |JoblibMinVersion|)
- threadpoolctl (>= |ThreadpoolctlMinVersion|)

Scikit-learn plotting capabilities (i.e., functions start with ``plot_`` and
classes end with ``Display``) require Matplotlib (>= |MatplotlibMinVersion|).
For running the examples Matplotlib >= |MatplotlibMinVersion| is required.
A few examples require scikit-image >= |Scikit-ImageMinVersion|, a few examples
require pandas >= |PandasMinVersion|, some examples require seaborn >=
|SeabornMinVersion| and Plotly >= |PlotlyMinVersion|.

User installation
~~~~~~~~~~~~~~~~~

If you already have a working installation of NumPy and SciPy,
the easiest way to install scikit-learn is using ``pip``::

    pip install -U scikit-learn

or ``conda``::

    conda install -c conda-forge scikit-learn

The documentation includes more detailed `installation instructions <https://scikit-learn.org/stable/install.html>`_.

Changelog
---------

See the `changelog <https://scikit-learn.org/dev/whats_new.html>`__
for a history of notable changes to scikit-learn.

Development
-----------

We welcome new contributors of all experience levels. The scikit-learn
community goals are to be helpful, welcoming, and effective. The
`Development Guide <https://scikit-learn.org/stable/developers/index.html>`_
has detailed information about contributing code, documentation, tests, and
more. We've included some basic information in this README.

Important links
~~~~~~~~~~~~~~~

- Official source code repo: https://github.com/scikit-learn/scikit-learn
- Download releases: https://pypi.org/project/scikit-learn/
- Issue tracker: https://github.com/scikit-learn/scikit-learn/issues

Source code
~~~~~~~~~~~

You can check the latest sources with the command::

    git clone https://github.com/scikit-learn/scikit-learn.git

Contributing
~~~~~~~~~~~~

To learn more about making a contribution to scikit-learn, please see our
`Contributing guide
<https://scikit-learn.org/dev/developers/contributing.html>`_.

Testing
~~~~~~~

After installation, you can launch the test suite from outside the source
directory (you will need to have ``pytest`` >= |PytestMinVersion| installed)::

    pytest sklearn

See the web page
https://scikit-learn.org/dev/developers/contributing.html#testing-and-improving-test-coverage
for more information.

Random number generation can be controlled during testing by setting
the ``SKLEARN_SEED`` environment variable.

Submitting a Pull Request
~~~~~~~~~~~~~~~~~~~~~~~~~

Before opening a Pull Request, have a look at the
full Contributing page to make sure your code complies
with our guidelines:
https://scikit-learn.org/stable/developers/index.html

Project History
---------------

The project was started in 2007 by David Cournapeau as a Google Summer
of Code project, and since then many volunteers have contributed. See
the `About us <https://scikit-learn.org/dev/about.html#authors>`__ page
for a list of core contributors.

The project is currently maintained by a team of volunteers.

**Note**: ``scikit-learn`` was previously referred to as ``scikits.learn``.

Help and Support
----------------

Documentation
~~~~~~~~~~~~~

- HTML documentation (stable release): https://scikit-learn.org
- HTML documentation (development version): https://scikit-learn.org/dev/
- FAQ: https://scikit-learn.org/stable/faq.html

Communication
~~~~~~~~~~~~~

Main Channels
^^^^^^^^^^^^^

- **Website**: https://scikit-learn.org
- **Blog**: https://blog.scikit-learn.org
- **Mailing list**:
  https://mail.python.org/mailman/listinfo/scikit-learn

Developer & Support
^^^^^^^^^^^^^^^^^^^

- **GitHub Discussions**:
  https://github.com/scikit-learn/scikit-learn/discussions
- **Stack Overflow**:
  https://stackoverflow.com/questions/tagged/scikit-learn
- **Discord**: https://discord.gg/h9qyrK8Jc8

Social Media Platforms
^^^^^^^^^^^^^^^^^^^^^^

- **LinkedIn**:
  https://www.linkedin.com/company/scikit-learn
- **YouTube**:
  https://www.youtube.com/channel/UCJosFjYm0ZYVUARxuOZqnnw/playlists
- **Facebook**:
  https://www.facebook.com/scikitlearnofficial/
- **Instagram**:
  https://www.instagram.com/scikitlearnofficial/
- **TikTok**:
  https://www.tiktok.com/@scikit.learn
- **Bluesky**:
  https://bsky.app/profile/scikit-learn.org
- **Mastodon**:
  https://mastodon.social/@sklearn@fosstodon.org

Resources
^^^^^^^^^

- **Calendar**:
  https://blog.scikit-learn.org/calendar/
- **Logos & Branding**:
  https://github.com/scikit-learn/scikit-learn/tree/main/doc/logos

Citation
~~~~~~~~

If you use scikit-learn in a scientific publication, we would appreciate
citations:
https://scikit-learn.org/stable/about.html#citing-scikit-learn
