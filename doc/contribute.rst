===========
Development
===========

This project is a community effort, and everyone is welcomed to
contribute.

Developers web site
===================
This is the central web page for developers
http://sourceforge.net/apps/trac/scikit-learn/wiki

Code
====


Git repo
--------

You can check the sources with the command::
    
    git clone git://scikit-learn.git.sourceforge.net/gitroot/scikit-learn/scikit-learn

If you have contributed some code and would like to have write
privileges in subversion repository, please contact me (Fabian
Pedregosa <fabian.pedregosa@inria.fr>) and I'll give you write
privileges for the svn.


Patches
-------
Patches are the prefered way to contribute to a project if you do not
have write priviles.

Let's suppose that you have the latest sources for subversion and that
you just made some modifications that you'd like to share with the
world. You might proceed as:

1. Create a patch file. The command::

    git format-patch origin

will create a series of patch files with the changes you made with
the code base. 

2. Send that file to the mailing list or attach it to an
issue in the issue tracker and some devs will push that patch to the
main repository.

3. Wait for a reply. You should soon recive a reply on wether your
patch was committed.


EasyFix Issues
^^^^^^^^^^^^^^

The best way to get your feet wet is to pick up an issue from the
`issue tracker
<https://sourceforge.net/apps/trac/scikit-learn/report>`_ that are
labeled as EasyFix. This means that the knowledge needed to solve the
issue is low, but still you are helping the project and letting more
experienced developers concentrate on other issues.



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

