cimport cython

from .Visitor cimport (
    CythonTransform, VisitorTransform, TreeVisitor,
    ScopeTrackingTransform, EnvTransform)

# Don't include mixins, only the main classes.
#cdef class SkipDeclarations:

cdef class NormalizeTree(CythonTransform):
    cdef bint is_in_statlist
    cdef bint is_in_expr
    cpdef visit_StatNode(self, node, is_listcontainer=*)

cdef class PostParse(ScopeTrackingTransform):
    cdef dict specialattribute_handlers
    cdef size_t lambda_counter
    cdef size_t genexpr_counter
    cdef bint in_pattern_node
    cdef _visit_assignment_node(self, node, list expr_list)


#class PxdPostParse(CythonTransform, SkipDeclarations):
#class InterpretCompilerDirectives(CythonTransform, SkipDeclarations):
#class WithTransform(VisitorTransform, SkipDeclarations):
#class DecoratorTransform(CythonTransform, SkipDeclarations):

#class AnalyseDeclarationsTransform(EnvTransform):

cdef class AnalyseExpressionsTransform(CythonTransform):
    cdef list positions

cdef class ExpandInplaceOperators(EnvTransform):
    pass

cdef class AlignFunctionDefinitions(CythonTransform):
    cdef dict directives
    cdef set imported_names
    cdef object scope

@cython.final
cdef class YieldNodeCollector(TreeVisitor):
    cdef public list yields
    cdef public list returns
    cdef public list finallys
    cdef public list excepts
    cdef public bint has_return_value
    cdef public bint has_yield
    cdef public bint has_await
    cdef list excludes

@cython.final
cdef class MarkClosureVisitor(CythonTransform):
    cdef bint needs_closure
    cdef list excludes

@cython.final
cdef class CreateClosureClasses(CythonTransform):
    cdef list path
    cdef bint in_lambda
    cdef module_scope
    cdef generator_class

    cdef create_class_from_scope(self, node, target_module_scope, inner_node=*)
    cdef find_entries_used_in_closures(self, node)

#cdef class InjectGilHandling(VisitorTransform, SkipDeclarations):
#    cdef bint nogil

cdef class GilCheck(VisitorTransform):
    cdef list env_stack
    cdef int nogil_state
    cdef int nogil_state_at_current_gilstatnode
    cdef object in_lock_block

cdef class TransformBuiltinMethods(EnvTransform):
    cdef dict def_node_body_insertions
    cdef visit_cython_attribute(self, node)
