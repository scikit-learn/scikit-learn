"""
This module provides classes and functions for pass management.

There are two kinds of passes: transformations and analysis.
    * ModuleAnalysis, FunctionAnalysis and NodeAnalysis are to be
      subclassed by any pass that collects information about the AST.
    * gather is used to gather (!) the result of an analyses on an AST node.
    * Backend is to be sub-classed by any pass that dumps the AST content.
    * dump is used to dump (!) the AST using the given backend.
    * Transformation is to be sub-classed by any pass that updates the AST.
    * apply is used to apply (sic) a transformation on an AST node.
"""

import gast as ast
import os
import re
from collections.abc import Iterable
from time import time


def uncamel(name):
    """Transform CamelCase naming convention into C-ish convention."""
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()


class AnalysisContext:

    """
    Class that stores the hierarchy of node visited.

    Contains:
        * parent module
        * parent function
    """

    def __init__(self):
        self.module = None
        self.function = None


class ContextManagerMeta(type):
    def __getitem__(cls, Dependencies):
        class CustomContextManager(cls):
            if isinstance(Dependencies, Iterable):
                deps = Dependencies
            else:
                deps = Dependencies,
        return CustomContextManager

class ContextManager(metaclass=ContextManagerMeta):

    """
    Class to be inherited from to add automatic update of context.

    The optional analysis dependencies are listed in `dependencies'.
    """

    deps = ()

    def __init__(self):
        """ Create default context and save all dependencies. """
        self.verify_dependencies()

    def attach(self, pm, ctx=None):
        self.passmanager = pm
        self.ctx = ctx or AnalysisContext()

    def verify_dependencies(self):
        """
        Checks no analysis are called before a transformation,
        as the transformation could invalidate the analysis.
        """

        for i in range(1, len(self.deps)):
            assert(not (isinstance(self.deps[i], Transformation) and
                   isinstance(self.deps[i - 1], Analysis))
                   ), "invalid dep order for %s" % self

    def visit(self, node):
        if isinstance(node, ast.FunctionDef):
            self.ctx.function = node
            for D in self.deps:
                if issubclass(D, FunctionAnalysis):
                    # this should have already been computed as part of the run
                    # method of function analysis triggered by prepare
                    result = self.passmanager._cache[node, D]
                    setattr(self, uncamel(D.__name__), result)
        return super(ContextManager, self).visit(node)

    def prepare(self, node):
        '''Gather analysis result required by this analysis'''
        if isinstance(node, ast.Module):
            self.ctx.module = node
        elif isinstance(node, ast.FunctionDef):
            self.ctx.function = node

        for D in self.deps:
            d = D()
            d.attach(self.passmanager, self.ctx)
            result = d.run(node)
            setattr(self, uncamel(D.__name__), result)

    def run(self, node):
        """Override this to add special pre or post processing handlers."""
        self.prepare(node)
        return self.visit(node)

    def gather(self, analysis, node):
        a = analysis()
        a.attach(self.passmanager, self.ctx)
        return a.run(node)

class Analysis(ContextManager, ast.NodeVisitor):
    """
    A pass that does not change its content but gathers informations about it.
    """
    def __init__(self):
        super().__init__()
        self.result = type(self).ResultType()
        self.update = False

    def run(self, node):
        key = node, type(self)
        if key in self.passmanager._cache:
            self.result = self.passmanager._cache[key]
        else:
            super(Analysis, self).run(node)
            self.passmanager._cache[key] = self.result
        return self.result

    def display(self, data):
        print(data)

    def apply(self, node):
        self.display(self.run(node))
        return False, node

class ModuleAnalysis(Analysis):
    """An analysis that operates on a whole module."""

    def run(self, node):
        if not isinstance(node, ast.Module):
            if self.ctx.module is None:
                raise ValueError("{} called in an uninitialized context"
                                 .format(type(self).__name__))
            node = self.ctx.module
        return super(ModuleAnalysis, self).run(node)

class FunctionAnalysis(Analysis):
    """An analysis that operates on a function."""

    def run(self, node):
        if isinstance(node, ast.Module):
            self.ctx.module = node
            last = None
            for stmt in node.body:
                if isinstance(stmt, ast.FunctionDef):
                    last = self.gather(type(self), stmt)
            # last is None if there's no function to process
            return self.result if last is None else last
        elif not isinstance(node, ast.FunctionDef):
            if self.ctx.function is None:
                raise ValueError("{} called in an uninitialized context"
                                 .format(type(self).__name__))
            node = self.ctx.function
        return super(FunctionAnalysis, self).run(node)

class NodeAnalysis(Analysis):
    """An analysis that operates on any node."""

class Backend(ModuleAnalysis):
    """A pass that produces code from an AST."""


class Transformation(ContextManager, ast.NodeTransformer):

    """A pass that updates its content."""

    def __init__(self):
        """ Initialize the update used to know if update happened. """
        super(Transformation, self).__init__()
        self.update = False

    def run(self, node):
        """ Apply transformation and dependencies and fix new node location."""
        n = super(Transformation, self).run(node)
        # the transformation updated the AST, so analyse may need to be rerun
        # we could use a finer-grain caching system, and provide a way to flag
        # some analyses as `unmodified' by the transformation, as done in LLVM
        # (and PIPS ;-)
        if self.update:
            self.passmanager._cache.clear()
        return n

    def apply(self, node):
        """ Apply transformation and return if an update happened. """
        new_node = self.run(node)
        return self.update, new_node


class PassManager:
    '''
    Front end to the pythran pass system.
    '''
    def __init__(self, module_name, module_dir=None, code=None):
        self.module_name = module_name
        self.module_dir = module_dir or os.getcwd()
        self.code = code
        self._cache = {}

    def gather(self, analysis, node, run_times = None):
        "High-level function to call an `analysis' on a `node'"
        t0 = time()
        assert issubclass(analysis, Analysis)
        a = analysis()
        a.attach(self)
        ret = a.run(node)
        if run_times is not None: run_times[analysis] = run_times.get(analysis,0) + time()-t0

        return ret

    def dump(self, backend, node):
        '''High-level function to call a `backend' on a `node' to generate
        code for module `module_name'.'''
        assert issubclass(backend, Backend)
        b = backend()
        b.attach(self)
        return b.run(node)

    def apply(self, transformation, node, run_times = None):
        '''
        High-level function to call a `transformation' on a `node'.
        If the transformation is an analysis, the result of the analysis
        is displayed.
        '''
        t0 = time()
        assert issubclass(transformation, (Transformation, Analysis))
        a = transformation()
        a.attach(self)
        ret=a.apply(node)
        if run_times is not None: run_times[transformation] = run_times.get(transformation,0) + time()-t0
        return ret
