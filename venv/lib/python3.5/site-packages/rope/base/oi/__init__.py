"""Rope object analysis and inference package

Rope makes some simplifying assumptions about a python program.  It
assumes that a program only performs assignments and function calls.
Tracking assignments is simple and `PyName` objects handle that.  The
main problem is function calls.  Rope uses these two approaches for
obtaining call information:

* Static object analysis: `rope.base.pycore.PyCore.analyze_module()`

  It can analyze modules to obtain information about functions.  This
  is done by analyzing function calls in a module or scope.  Currently
  SOA analyzes the scopes that are changed while saving or when the
  user asks to analyze a module.  That is mainly because static
  analysis is time-consuming.

* Dynamic object analysis: `rope.base.pycore.PyCore.run_module()`

  When you run a module or your testsuite, when DOA is enabled, it
  collects information about parameters passed to and objects returned
  from functions.  The main problem with this approach is that it is
  quite slow; Not when looking up the information but when collecting
  them.

An instance of `rope.base.oi.objectinfo.ObjectInfoManager` can be used
for accessing these information.  It saves the data in a
`rope.base.oi.objectdb.ObjectDB` internally.

Now if our objectdb does not know anything about a function and we
need the value returned by it, static object inference, SOI, comes
into play.  It analyzes function body and tries to infer the object
that is returned from it (we usually need the returned value for the
given parameter objects).

Rope might collect and store information for other `PyName`\s, too.
For instance rope stores the object builtin containers hold.

"""
