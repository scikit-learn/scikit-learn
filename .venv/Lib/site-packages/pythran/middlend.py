"""This module turns a python AST into an optimized, pythran compatible ast."""

from pythran.analyses import ExtendedSyntaxCheck
from pythran.optimizations import (ComprehensionPatterns, ListCompToGenexp,
                                   RemoveDeadFunctions)
from pythran.transformations import (ExpandBuiltins, ExpandImports,
                                     ExpandImportAll, FalsePolymorphism,
                                     NormalizeCompare, NormalizeException,
                                     NormalizeMethodCalls, NormalizeReturn,
                                     NormalizeTuples, RemoveComprehension,
                                     RemoveNestedFunctions, RemoveLambdas,
                                     UnshadowParameters, RemoveNamedArguments,
                                     ExpandGlobals, NormalizeIsNone,
                                     NormalizeIfElse, NormalizeTypeIs,
                                     NormalizeStaticIf, SplitStaticExpression,
                                     RemoveFStrings)
import re


def refine(pm, node, optimizations, report_times=False):
    """ Refine node in place until it matches pythran's expectations. """
    run_times = {}
    # Sanitize input
    pm.apply(RemoveDeadFunctions, node, run_times)
    pm.apply(ExpandGlobals, node, run_times)
    pm.apply(ExpandImportAll, node, run_times)
    pm.apply(NormalizeTuples, node, run_times)
    pm.apply(RemoveFStrings, node, run_times)
    pm.apply(ExpandBuiltins, node, run_times)
    pm.apply(ExpandImports, node, run_times)
    pm.apply(NormalizeMethodCalls, node, run_times)
    pm.apply(NormalizeTypeIs, node, run_times)
    pm.apply(NormalizeIfElse, node, run_times)
    pm.apply(NormalizeIsNone, node, run_times)
    pm.apply(SplitStaticExpression, node, run_times)
    pm.apply(NormalizeStaticIf, node, run_times)
    pm.apply(NormalizeTuples, node, run_times)
    pm.apply(NormalizeException, node, run_times)
    pm.apply(NormalizeMethodCalls, node, run_times)

    # Some early optimizations
    pm.apply(ComprehensionPatterns, node, run_times)

    pm.apply(RemoveLambdas, node, run_times)
    pm.apply(RemoveNestedFunctions, node, run_times)
    pm.apply(NormalizeCompare, node, run_times)

    pm.gather(ExtendedSyntaxCheck, node, run_times)

    pm.apply(ListCompToGenexp, node, run_times)
    pm.apply(RemoveComprehension, node, run_times)
    pm.apply(RemoveNamedArguments, node, run_times)

    # sanitize input
    pm.apply(NormalizeReturn, node, run_times)
    pm.apply(UnshadowParameters, node, run_times)
    pm.apply(FalsePolymorphism, node, run_times)

    # some extra optimizations
    apply_optimisation = True
    while apply_optimisation:
        apply_optimisation = False
        for optimization in optimizations:
            apply_optimisation |= pm.apply(optimization, node, run_times)[0]

    if report_times and len(run_times):
        print("Optimization run times:")
        for key,val in run_times.items():
            k = re.findall("'(.*)'",str(key))
            print(f"{k[0]:<70s} : {val:>5.1f} s")

def mark_unexported_functions(ir, exported_functions):
    from pythran.metadata import add as MDadd, Local as MDLocal
    for node in ir.body:
        if hasattr(node, 'name'):
            if node.name not in exported_functions:
                MDadd(node, MDLocal())
    return ir
