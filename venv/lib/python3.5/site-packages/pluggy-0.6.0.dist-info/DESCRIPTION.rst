pluggy - A minimalist production ready plugin system
====================================================
|pypi| |anaconda| |versions| |travis| |appveyor|


This is the core framework used by the `pytest`_, `tox`_, and `devpi`_ projects.

Please `read the docs`_ to learn more!

A definitive example
********************
.. code-block:: python

    import pluggy

    hookspec = pluggy.HookspecMarker("myproject")
    hookimpl = pluggy.HookimplMarker("myproject")


    class MySpec(object):
        """A hook specification namespace.
        """
        @hookspec
        def myhook(self, arg1, arg2):
            """My special little hook that you can customize.
            """


    class Plugin_1(object):
        """A hook implementation namespace.
        """
        @hookimpl
        def myhook(self, arg1, arg2):
            print("inside Plugin_1.myhook()")
            return arg1 + arg2


    class Plugin_2(object):
        """A 2nd hook implementation namespace.
        """
        @hookimpl
        def myhook(self, arg1, arg2):
            print("inside Plugin_2.myhook()")
            return arg1 - arg2


    # create a manager and add the spec
    pm = pluggy.PluginManager("myproject")
    pm.add_hookspecs(MySpec)

    # register plugins
    pm.register(Plugin_1())
    pm.register(Plugin_2())

    # call our `myhook` hook
    results = pm.hook.myhook(arg1=1, arg2=2)
    print(results)


.. badges
.. |pypi| image:: https://img.shields.io/pypi/v/pluggy.svg
    :target: https://pypi.python.org/pypi/pluggy
.. |versions| image:: https://img.shields.io/pypi/pyversions/pluggy.svg
    :target: https://pypi.python.org/pypi/pluggy
.. |travis| image:: https://img.shields.io/travis/pytest-dev/pluggy/master.svg
    :target: https://travis-ci.org/pytest-dev/pluggy
.. |appveyor| image:: https://img.shields.io/appveyor/ci/pytestbot/pluggy/master.svg
    :target: https://ci.appveyor.com/project/pytestbot/pluggy
.. |anaconda| image:: https://anaconda.org/conda-forge/pluggy/badges/version.svg
    :target: https://anaconda.org/conda-forge/pluggy

.. links
.. _pytest:
    http://pytest.org
.. _tox:
    https://tox.readthedocs.org
.. _devpi:
    http://doc.devpi.net
.. _read the docs:
   https://pluggy.readthedocs.io/en/latest/


