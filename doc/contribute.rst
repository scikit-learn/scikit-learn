Contribute
==========

This project is a community effort, so everyone is welcomed to
contribute.


Code
----

EasyFix Issues
^^^^^^^^^^^^^^

The best way to get your feet wet is to pick up an issue from the
`issue tracker
<https://sourceforge.net/apps/trac/scikit-learn/report>`_ that are
labeled as EasyFix. This means that the knowledge needed to solve the
issue is low, but still you are helping the project and letting more
experienced developers concentrate on other issues.


Other Issues
^^^^^^^^^^^^

Just pick up an issue that is not assigned from the `issue tracker
<https://sourceforge.net/apps/trac/scikit-learn/report>`_ and submit a patch.



Patches
-------

Suppose you have made some modifications to the source code and want
to share your changes. The best way to share those changes is by
submitting a patch file. I suppose that you have checked the latest
sources from subversion. If you run the command::

  svn diff > patch.diff

This will create a file "patch.diff" with the changes you made with
the code base. Send that file to the mailing list or attach it to an
issue in the issue tracker and some devs will push that patch to the
main repository.

For more info about Subversion, you can read the excellent book
`Version Control with Subversion <http://svnbook.red-bean.com/>`_


SVN Access
----------

If you have contributed some code and would like to have write
privileges in subversion repository, please contact me (Fabian
Pedregosa <fabian.pedregosa@inria.fr>) and I'll add you to the list.


Git repo
--------

Some people find easier to work with Decentralized Version Control
Systems like Git. If that is your case, you can use a Git mirror that
is usually up to date with the main subversion repo. It's web
interface is located `here <http://github.com/fseoane/scikit-learn>`_
and the clone command is::

  git clone http://github.com/fseoane/scikit-learn


.. _packaging:

Packaging
^^^^^^^^^

You can also help making binary distributions for windows, OsX or packages for some
distribution.


Documentation
-------------

You can also contribute with documentation: function docstrings, rst
docs (like this one), tutorials, etc. Rst docs live in the source code
repository, under directory doc/.

