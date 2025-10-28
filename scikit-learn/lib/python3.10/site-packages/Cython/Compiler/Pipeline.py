import itertools
from time import time

from . import Errors
from . import DebugFlags
from . import Options
from .Errors import CompileError, InternalError, AbortError
from . import Naming

#
# Really small pipeline stages
#
def dumptree(t):
    # For quick debugging in pipelines
    print(t.dump())
    return t

def abort_on_errors(node):
    # Stop the pipeline if there are any errors.
    if Errors.get_errors_count() != 0:
        raise AbortError("pipeline break")
    return node

def parse_stage_factory(context):
    def parse(compsrc):
        source_desc = compsrc.source_desc
        full_module_name = compsrc.full_module_name
        initial_pos = (source_desc, 1, 0)
        saved_cimport_from_pyx, Options.cimport_from_pyx = Options.cimport_from_pyx, False
        scope = context.find_module(full_module_name, pos = initial_pos, need_pxd = 0)
        Options.cimport_from_pyx = saved_cimport_from_pyx
        tree = context.parse(source_desc, scope, pxd = 0, full_module_name = full_module_name)
        tree.compilation_source = compsrc
        tree.scope = scope
        tree.is_pxd = False
        return tree
    return parse

def parse_pxd_stage_factory(context, scope, module_name):
    def parse(source_desc):
        tree = context.parse(source_desc, scope, pxd=True,
                             full_module_name=module_name)
        tree.scope = scope
        tree.is_pxd = True
        return tree
    return parse

def generate_pyx_code_stage_factory(options, result):
    def generate_pyx_code_stage(module_node):
        module_node.process_implementation(options, result)
        result.compilation_source = module_node.compilation_source
        return result
    return generate_pyx_code_stage

def inject_utility_pxd_code_stage_factory(context):

    def inject_utility_pxd_code_stage(module_node):
        for statlistnode, scope in context.utility_pxds.values():
            module_node.merge_in(statlistnode, scope, stage="pxd")
        return module_node

    return inject_utility_pxd_code_stage

def inject_pxd_code_stage_factory(context):

    def inject_pxd_code_stage(module_node):
        for name, (statlistnode, scope) in context.pxds.items():
            module_node.merge_in(statlistnode, scope, stage="pxd")
        return module_node
    return inject_pxd_code_stage


def use_utility_code_definitions(scope, target, seen=None):
    if seen is None:
        seen = set()

    for entry in scope.entries.values():
        if entry in seen:
            continue

        seen.add(entry)
        if entry.used and entry.utility_code_definition:
            target.use_utility_code(entry.utility_code_definition)
            for required_utility in entry.utility_code_definition.requires:
                target.use_utility_code(required_utility)
        elif entry.as_module:
            use_utility_code_definitions(entry.as_module, target, seen)


def sorted_utility_codes_and_deps(utilcodes):
    ranks = {}
    get_rank = ranks.get

    def calculate_rank(utilcode):
        rank = get_rank(utilcode)
        if rank is None:
            ranks[utilcode] = 0  # prevent infinite recursion on circular dependencies
            original_order = len(ranks)
            rank = ranks[utilcode] = 1 + (
                min([calculate_rank(dep) for dep in utilcode.requires]) if utilcode.requires else -1
                ) + original_order * 1e-8
        return rank

    for utilcode in utilcodes:
        calculate_rank(utilcode)

    # include all recursively collected dependencies
    return sorted(ranks, key=get_rank)


def normalize_deps(utilcodes):
    deps = {utilcode:utilcode for utilcode in utilcodes}
    for utilcode in utilcodes:
        utilcode.requires = [deps.setdefault(dep, dep) for dep in utilcode.requires or ()]


def inject_utility_code_stage_factory(context, internalise_c_class_entries=True):
    def inject_utility_code_stage(module_node):
        module_node.prepare_utility_code()
        use_utility_code_definitions(context.cython_scope, module_node.scope)

        module_scope = module_node.scope
        utility_code_list = module_scope.utility_code_list
        utility_code_list[:] = sorted_utility_codes_and_deps(utility_code_list)
        normalize_deps(utility_code_list)

        added = set()
        # Note: the list might be extended inside the loop (if some utility code
        # pulls in other utility code, explicitly or implicitly)
        for utilcode in utility_code_list:
            if utilcode in added:
                continue
            added.add(utilcode)
            if utilcode.requires:
                for dep in utilcode.requires:
                    if dep not in added:
                        utility_code_list.append(dep)
            if tree := utilcode.get_tree(cython_scope=context.cython_scope):
                module_node.merge_in(tree.with_compiler_directives(),
                                     tree.scope, stage="utility")
                module_node.merge_scope(tree.scope, internalise_c_class_entries=internalise_c_class_entries)
            elif shared_library_scope := utilcode.get_shared_library_scope(cython_scope=context.cython_scope):
                module_scope.cimported_modules.append(shared_library_scope)
        return module_node

    return inject_utility_code_stage


#
# Pipeline factories
#

def create_pipeline(context, mode, exclude_classes=()):
    assert mode in ('pyx', 'py', 'pxd')
    from .Visitor import PrintTree
    from .ParseTreeTransforms import WithTransform, NormalizeTree, PostParse, PxdPostParse
    from .ParseTreeTransforms import ForwardDeclareTypes, InjectGilHandling, AnalyseDeclarationsTransform
    from .ParseTreeTransforms import AnalyseExpressionsTransform, FindInvalidUseOfFusedTypes
    from .ParseTreeTransforms import CreateClosureClasses, MarkClosureVisitor, DecoratorTransform
    from .ParseTreeTransforms import TrackNumpyAttributes, InterpretCompilerDirectives, TransformBuiltinMethods
    from .ParseTreeTransforms import ExpandInplaceOperators, ParallelRangeTransform
    from .ParseTreeTransforms import CalculateQualifiedNamesTransform
    from .TypeInference import MarkParallelAssignments, MarkOverflowingArithmetic
    from .ParseTreeTransforms import AdjustDefByDirectives, AlignFunctionDefinitions, AutoCpdefFunctionDefinitions
    from .ParseTreeTransforms import RemoveUnreachableCode, GilCheck, CoerceCppTemps
    from .FlowControl import ControlFlowAnalysis
    from .AnalysedTreeTransforms import AutoTestDictTransform
    from .AutoDocTransforms import EmbedSignature
    from .Optimize import FlattenInListTransform, SwitchTransform, IterationTransform
    from .Optimize import EarlyReplaceBuiltinCalls, OptimizeBuiltinCalls
    from .Optimize import InlineDefNodeCalls
    from .Optimize import ConstantFolding, FinalOptimizePhase
    from .Optimize import DropRefcountingTransform
    from .Optimize import ConsolidateOverflowCheck
    from .Buffer import IntroduceBufferAuxiliaryVars
    from .ModuleNode import check_c_declarations, check_c_declarations_pxd


    if mode == 'pxd':
        _check_c_declarations = check_c_declarations_pxd
        _specific_post_parse = PxdPostParse(context)
    else:
        _check_c_declarations = check_c_declarations
        _specific_post_parse = None

    if mode == 'py':
        _align_function_definitions = AlignFunctionDefinitions(context)
    else:
        _align_function_definitions = None

    # NOTE: This is the "common" parts of the pipeline, which is also
    # code in pxd files. So it will be run multiple times in a
    # compilation stage.
    stages = [
        NormalizeTree(context),
        PostParse(context),
        _specific_post_parse,
        TrackNumpyAttributes(),
        InterpretCompilerDirectives(context, context.compiler_directives),
        ParallelRangeTransform(context),
        WithTransform(),
        AdjustDefByDirectives(context),
        _align_function_definitions,
        MarkClosureVisitor(context),
        AutoCpdefFunctionDefinitions(context),
        RemoveUnreachableCode(context),
        ConstantFolding(),
        FlattenInListTransform(),
        DecoratorTransform(context),
        ForwardDeclareTypes(context),
        InjectGilHandling(),
        AnalyseDeclarationsTransform(context),
        AutoTestDictTransform(context),
        EmbedSignature(context),
        EarlyReplaceBuiltinCalls(context),  ## Necessary?
        TransformBuiltinMethods(context),
        MarkParallelAssignments(context),
        ControlFlowAnalysis(context),
        RemoveUnreachableCode(context),
        # MarkParallelAssignments(context),
        MarkOverflowingArithmetic(context),
        IntroduceBufferAuxiliaryVars(context),
        _check_c_declarations,
        InlineDefNodeCalls(context),
        AnalyseExpressionsTransform(context),
        FindInvalidUseOfFusedTypes(),
        ExpandInplaceOperators(context),
        IterationTransform(context),
        SwitchTransform(context),
        OptimizeBuiltinCalls(context),  ## Necessary?
        CreateClosureClasses(context),  ## After all lookups and type inference
        CalculateQualifiedNamesTransform(context),
        ConsolidateOverflowCheck(context),
        DropRefcountingTransform(),
        FinalOptimizePhase(context),
        CoerceCppTemps(context),
        GilCheck(),
        ]
    if exclude_classes:
        stages = [s for s in stages if s.__class__ not in exclude_classes]
    return stages

def create_pyx_pipeline(context, options, result, py=False, exclude_classes=()):
    mode = 'py' if py else 'pyx'

    test_support = []
    ctest_support = []
    if options.evaluate_tree_assertions:
        from ..TestUtils import TreeAssertVisitor
        test_validator = TreeAssertVisitor()
        test_support.append(test_validator)
        ctest_support.append(test_validator.create_c_file_validator())

    if options.gdb_debug:
        from ..Debugger import DebugWriter  # requires Py2.5+
        from .ParseTreeTransforms import DebugTransform
        context.gdb_debug_outputwriter = DebugWriter.CythonDebugWriter(
            options.output_dir)
        debug_transform = [DebugTransform(context, options, result)]
    else:
        debug_transform = []

    return list(itertools.chain(
        [parse_stage_factory(context)],
        create_pipeline(context, mode, exclude_classes=exclude_classes),
        test_support,
        [
            inject_pxd_code_stage_factory(context),
            inject_utility_code_stage_factory(context),
            inject_utility_pxd_code_stage_factory(context),
            abort_on_errors,
        ],
        debug_transform,
        [generate_pyx_code_stage_factory(options, result)],
        ctest_support,
    ))

def create_pxd_pipeline(context, scope, module_name):
    from .CodeGeneration import ExtractPxdCode

    # The pxd pipeline ends up with a CCodeWriter containing the
    # code of the pxd, as well as a pxd scope.
    return [
        parse_pxd_stage_factory(context, scope, module_name)
        ] + create_pipeline(context, 'pxd') + [
        ExtractPxdCode()
        ]

def create_py_pipeline(context, options, result):
    return create_pyx_pipeline(context, options, result, py=True)

def create_pyx_as_pxd_pipeline(context, result):
    from .ParseTreeTransforms import AlignFunctionDefinitions, \
        MarkClosureVisitor, WithTransform, AnalyseDeclarationsTransform
    from .Optimize import ConstantFolding, FlattenInListTransform
    from .Nodes import StatListNode
    pipeline = []
    pyx_pipeline = create_pyx_pipeline(context, context.options, result,
                                       exclude_classes=[
                                           AlignFunctionDefinitions,
                                           MarkClosureVisitor,
                                           ConstantFolding,
                                           FlattenInListTransform,
                                           WithTransform
                                           ])
    from .Visitor import VisitorTransform
    class SetInPxdTransform(VisitorTransform):
        # A number of nodes have an "in_pxd" attribute which affects AnalyseDeclarationsTransform
        # (for example controlling pickling generation). Set it, to make sure we don't mix them up with
        # the importing main module.
        # FIXME: This should be done closer to the parsing step.
        def visit_StatNode(self, node):
            if hasattr(node, "in_pxd"):
                node.in_pxd = True
            self.visitchildren(node)
            return node

        visit_Node = VisitorTransform.recurse_to_children

    for stage in pyx_pipeline:
        pipeline.append(stage)
        if isinstance(stage, AnalyseDeclarationsTransform):
            pipeline.insert(-1, SetInPxdTransform())
            break  # This is the last stage we need.
    def fake_pxd(root):
        for entry in root.scope.entries.values():
            if not entry.in_cinclude:
                entry.defined_in_pxd = 1
                if entry.name == entry.cname and entry.visibility != 'extern':
                    # Always mangle non-extern cimported entries.
                    entry.cname = entry.scope.mangle(Naming.func_prefix, entry.name)
        return StatListNode(root.pos, stats=[]), root.scope
    pipeline.append(fake_pxd)
    return pipeline

def insert_into_pipeline(pipeline, transform, before=None, after=None):
    """
    Insert a new transform into the pipeline after or before an instance of
    the given class. e.g.

        pipeline = insert_into_pipeline(pipeline, transform,
                                        after=AnalyseDeclarationsTransform)
    """
    assert before or after

    cls = before or after
    for i, t in enumerate(pipeline):
        if isinstance(t, cls):
            break

    if after:
        i += 1

    return pipeline[:i] + [transform] + pipeline[i:]

#
# Running a pipeline
#

try:
    from threading import local as _threadlocal
except ImportError:
    class _threadlocal: pass

threadlocal = _threadlocal()


def get_timings():
    try:
        return threadlocal.cython_pipeline_timings
    except AttributeError:
        return {}


_pipeline_entry_points = {}

def _make_debug_phase_runner(phase_name):
    # Create a new wrapper for each step to show the name in profiles.
    try:
        return _pipeline_entry_points[phase_name]
    except KeyError:
        pass

    def run(phase, data):
        return phase(data)

    run.__name__ = run.__qualname__ = phase_name
    _pipeline_entry_points[phase_name] = run
    return run


def run_pipeline(pipeline, source, printtree=True):
    from .Visitor import PrintTree
    try:
        timings = threadlocal.cython_pipeline_timings
    except AttributeError:
        timings = threadlocal.cython_pipeline_timings = {}

    def run(phase, data):
        return phase(data)

    error = None
    data = source
    try:
        try:
            for phase in pipeline:
                if phase is None:
                    continue
                if not printtree and isinstance(phase, PrintTree):
                    continue

                phase_name = getattr(phase, '__name__', type(phase).__name__)
                if DebugFlags.debug_verbose_pipeline:
                    print("Entering pipeline phase %r" % phase)
                    run = _make_debug_phase_runner(phase_name)

                t = time()
                data = run(phase, data)
                t = time() - t

                try:
                    old_t, count = timings[phase_name]
                except KeyError:
                    old_t, count = 0, 0
                timings[phase_name] = (old_t + int(t * 1000000), count + 1)
                if DebugFlags.debug_verbose_pipeline:
                    print("    %.3f seconds" % t)
        except CompileError as err:
            # err is set
            Errors.report_error(err, use_stack=False)
            error = err
    except InternalError as err:
        # Only raise if there was not an earlier error
        if Errors.get_errors_count() == 0:
            raise
        error = err
    except AbortError as err:
        error = err
    return (error, data)
