==========
Contribute
==========

This project is a community effort, and everyone is welcomed to
contribute.


Code
====

Patches
-------
Patches are the prefered way to contribute to a project if you do not
have (yet) write priviles.

Let's suppose that you have the latest sources for subversion and that
you just made some modifications that you'd like to share with the
world. The way to proceed is the following:

1. Create a patch file. The command::
    svn diff > patch.diff

will create a file "patch.diff" with the changes you made with
the code base. 

2. Send that file to the mailing list or attach it to an
issue in the issue tracker and some devs will push that patch to the
main repository.

3. Wait for a reply. You should soon recive a reply on wether your
patch was committed.

For more info about Subversion, you can read the excellent book
`Version Control with Subversion <http://svnbook.red-bean.com/>`_


EasyFix Issues
^^^^^^^^^^^^^^

The best way to get your feet wet is to pick up an issue from the
`issue tracker
<https://sourceforge.net/apps/trac/scikit-learn/report>`_ that are
labeled as EasyFix. This means that the knowledge needed to solve the
issue is low, but still you are helping the project and letting more
experienced developers concentrate on other issues.



SVN Access
----------

If you have contributed some code and would like to have write
privileges in subversion repository, please contact me (Fabian
Pedregosa <fabian.pedregosa@inria.fr>) and I'll give you write
privileges for the svn.


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
=============

I am glad to accept any sort of documentation: function docstrings, rst docs (like
this one), tutorials, etc. Rst docs live in the source code
repository, under directory doc/.

