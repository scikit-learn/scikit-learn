"""
Evaluation of Python code in |jedi| is based on three assumptions:

* The code uses as least side effects as possible. Jedi understands certain
  list/tuple/set modifications, but there's no guarantee that Jedi detects
  everything (list.append in different modules for example).
* No magic is being used:

  - metaclasses
  - ``setattr()`` / ``__import__()``
  - writing to ``globals()``, ``locals()``, ``object.__dict__``
* The programmer is not a total dick, e.g. like `this
  <https://github.com/davidhalter/jedi/issues/24>`_ :-)

The actual algorithm is based on a principle called lazy evaluation.  That
said, the typical entry point for static analysis is calling
``eval_expr_stmt``. There's separate logic for autocompletion in the API, the
evaluator is all about evaluating an expression.

TODO this paragraph is not what jedi does anymore, it's similar, but not the
same.

Now you need to understand what follows after ``eval_expr_stmt``. Let's
make an example::

    import datetime
    datetime.date.toda# <-- cursor here

First of all, this module doesn't care about completion. It really just cares
about ``datetime.date``. At the end of the procedure ``eval_expr_stmt`` will
return the ``date`` class.

To *visualize* this (simplified):

- ``Evaluator.eval_expr_stmt`` doesn't do much, because there's no assignment.
- ``Context.eval_node`` cares for resolving the dotted path
- ``Evaluator.find_types`` searches for global definitions of datetime, which
  it finds in the definition of an import, by scanning the syntax tree.
- Using the import logic, the datetime module is found.
- Now ``find_types`` is called again by ``eval_node`` to find ``date``
  inside the datetime module.

Now what would happen if we wanted ``datetime.date.foo.bar``? Two more
calls to ``find_types``. However the second call would be ignored, because the
first one would return nothing (there's no foo attribute in ``date``).

What if the import would contain another ``ExprStmt`` like this::

    from foo import bar
    Date = bar.baz

Well... You get it. Just another ``eval_expr_stmt`` recursion. It's really
easy. Python can obviously get way more complicated then this. To understand
tuple assignments, list comprehensions and everything else, a lot more code had
to be written.

Jedi has been tested very well, so you can just start modifying code. It's best
to write your own test first for your "new" feature. Don't be scared of
breaking stuff. As long as the tests pass, you're most likely to be fine.

I need to mention now that lazy evaluation is really good because it
only *evaluates* what needs to be *evaluated*. All the statements and modules
that are not used are just being ignored.
"""

from parso.python import tree
import parso
from parso import python_bytes_to_unicode

from jedi import debug
from jedi import parser_utils
from jedi.evaluate.utils import unite
from jedi.evaluate import imports
from jedi.evaluate import recursion
from jedi.evaluate.cache import evaluator_function_cache
from jedi.evaluate import compiled
from jedi.evaluate import helpers
from jedi.evaluate.filters import TreeNameDefinition, ParamName
from jedi.evaluate.base_context import ContextualizedName, ContextualizedNode, \
    ContextSet, NO_CONTEXTS, iterate_contexts
from jedi.evaluate.context import ClassContext, FunctionContext, \
    AnonymousInstance, BoundMethod
from jedi.evaluate.context.iterable import CompForContext
from jedi.evaluate.syntax_tree import eval_trailer, eval_expr_stmt, \
    eval_node, check_tuple_assignments


class Evaluator(object):
    def __init__(self, project, environment=None, script_path=None):
        if environment is None:
            environment = project.get_environment()
        self.environment = environment
        self.script_path = script_path
        self.compiled_subprocess = environment.get_evaluator_subprocess(self)
        self.grammar = environment.get_grammar()

        self.latest_grammar = parso.load_grammar(version='3.6')
        self.memoize_cache = {}  # for memoize decorators
        self.module_cache = imports.ModuleCache()  # does the job of `sys.modules`.
        self.compiled_cache = {}  # see `evaluate.compiled.create()`
        self.inferred_element_counts = {}
        self.mixed_cache = {}  # see `evaluate.compiled.mixed._create()`
        self.analysis = []
        self.dynamic_params_depth = 0
        self.is_analysis = False
        self.project = project
        self.access_cache = {}

        self.reset_recursion_limitations()
        self.allow_different_encoding = True

    @property
    @evaluator_function_cache()
    def builtins_module(self):
        return compiled.get_special_object(self, u'BUILTINS')

    def reset_recursion_limitations(self):
        self.recursion_detector = recursion.RecursionDetector()
        self.execution_recursion_detector = recursion.ExecutionRecursionDetector(self)

    def get_sys_path(self):
        """Convenience function"""
        return self.project._get_sys_path(self, environment=self.environment)

    def eval_element(self, context, element):
        if isinstance(context, CompForContext):
            return eval_node(context, element)

        if_stmt = element
        while if_stmt is not None:
            if_stmt = if_stmt.parent
            if if_stmt.type in ('if_stmt', 'for_stmt'):
                break
            if parser_utils.is_scope(if_stmt):
                if_stmt = None
                break
        predefined_if_name_dict = context.predefined_names.get(if_stmt)
        # TODO there's a lot of issues with this one. We actually should do
        # this in a different way. Caching should only be active in certain
        # cases and this all sucks.
        if predefined_if_name_dict is None and if_stmt \
                and if_stmt.type == 'if_stmt' and self.is_analysis:
            if_stmt_test = if_stmt.children[1]
            name_dicts = [{}]
            # If we already did a check, we don't want to do it again -> If
            # context.predefined_names is filled, we stop.
            # We don't want to check the if stmt itself, it's just about
            # the content.
            if element.start_pos > if_stmt_test.end_pos:
                # Now we need to check if the names in the if_stmt match the
                # names in the suite.
                if_names = helpers.get_names_of_node(if_stmt_test)
                element_names = helpers.get_names_of_node(element)
                str_element_names = [e.value for e in element_names]
                if any(i.value in str_element_names for i in if_names):
                    for if_name in if_names:
                        definitions = self.goto_definitions(context, if_name)
                        # Every name that has multiple different definitions
                        # causes the complexity to rise. The complexity should
                        # never fall below 1.
                        if len(definitions) > 1:
                            if len(name_dicts) * len(definitions) > 16:
                                debug.dbg('Too many options for if branch evaluation %s.', if_stmt)
                                # There's only a certain amount of branches
                                # Jedi can evaluate, otherwise it will take to
                                # long.
                                name_dicts = [{}]
                                break

                            original_name_dicts = list(name_dicts)
                            name_dicts = []
                            for definition in definitions:
                                new_name_dicts = list(original_name_dicts)
                                for i, name_dict in enumerate(new_name_dicts):
                                    new_name_dicts[i] = name_dict.copy()
                                    new_name_dicts[i][if_name.value] = ContextSet(definition)

                                name_dicts += new_name_dicts
                        else:
                            for name_dict in name_dicts:
                                name_dict[if_name.value] = definitions
            if len(name_dicts) > 1:
                result = ContextSet()
                for name_dict in name_dicts:
                    with helpers.predefine_names(context, if_stmt, name_dict):
                        result |= eval_node(context, element)
                return result
            else:
                return self._eval_element_if_evaluated(context, element)
        else:
            if predefined_if_name_dict:
                return eval_node(context, element)
            else:
                return self._eval_element_if_evaluated(context, element)

    def _eval_element_if_evaluated(self, context, element):
        """
        TODO This function is temporary: Merge with eval_element.
        """
        parent = element
        while parent is not None:
            parent = parent.parent
            predefined_if_name_dict = context.predefined_names.get(parent)
            if predefined_if_name_dict is not None:
                return eval_node(context, element)
        return self._eval_element_cached(context, element)

    @evaluator_function_cache(default=NO_CONTEXTS)
    def _eval_element_cached(self, context, element):
        return eval_node(context, element)

    def goto_definitions(self, context, name):
        def_ = name.get_definition(import_name_always=True)
        if def_ is not None:
            type_ = def_.type
            if type_ == 'classdef':
                return [ClassContext(self, context, name.parent)]
            elif type_ == 'funcdef':
                return [FunctionContext(self, context, name.parent)]

            if type_ == 'expr_stmt':
                is_simple_name = name.parent.type not in ('power', 'trailer')
                if is_simple_name:
                    return eval_expr_stmt(context, def_, name)
            if type_ == 'for_stmt':
                container_types = context.eval_node(def_.children[3])
                cn = ContextualizedNode(context, def_.children[3])
                for_types = iterate_contexts(container_types, cn)
                c_node = ContextualizedName(context, name)
                return check_tuple_assignments(self, c_node, for_types)
            if type_ in ('import_from', 'import_name'):
                return imports.infer_import(context, name)

        return helpers.evaluate_call_of_leaf(context, name)

    def goto(self, context, name):
        definition = name.get_definition(import_name_always=True)
        if definition is not None:
            type_ = definition.type
            if type_ == 'expr_stmt':
                # Only take the parent, because if it's more complicated than just
                # a name it's something you can "goto" again.
                is_simple_name = name.parent.type not in ('power', 'trailer')
                if is_simple_name:
                    return [TreeNameDefinition(context, name)]
            elif type_ == 'param':
                return [ParamName(context, name)]
            elif type_ in ('funcdef', 'classdef'):
                return [TreeNameDefinition(context, name)]
            elif type_ in ('import_from', 'import_name'):
                module_names = imports.infer_import(context, name, is_goto=True)
                return module_names

        par = name.parent
        node_type = par.type
        if node_type == 'argument' and par.children[1] == '=' and par.children[0] == name:
            # Named param goto.
            trailer = par.parent
            if trailer.type == 'arglist':
                trailer = trailer.parent
            if trailer.type != 'classdef':
                if trailer.type == 'decorator':
                    context_set = context.eval_node(trailer.children[1])
                else:
                    i = trailer.parent.children.index(trailer)
                    to_evaluate = trailer.parent.children[:i]
                    if to_evaluate[0] == 'await':
                        to_evaluate.pop(0)
                    context_set = context.eval_node(to_evaluate[0])
                    for trailer in to_evaluate[1:]:
                        context_set = eval_trailer(context, context_set, trailer)
                param_names = []
                for context in context_set:
                    try:
                        get_param_names = context.get_param_names
                    except AttributeError:
                        pass
                    else:
                        for param_name in get_param_names():
                            if param_name.string_name == name.value:
                                param_names.append(param_name)
                return param_names
        elif node_type == 'dotted_name':  # Is a decorator.
            index = par.children.index(name)
            if index > 0:
                new_dotted = helpers.deep_ast_copy(par)
                new_dotted.children[index - 1:] = []
                values = context.eval_node(new_dotted)
                return unite(
                    value.py__getattribute__(name, name_context=context, is_goto=True)
                    for value in values
                )

        if node_type == 'trailer' and par.children[0] == '.':
            values = helpers.evaluate_call_of_leaf(context, name, cut_own_trailer=True)
            return unite(
                value.py__getattribute__(name, name_context=context, is_goto=True)
                for value in values
            )
        else:
            stmt = tree.search_ancestor(
                name, 'expr_stmt', 'lambdef'
            ) or name
            if stmt.type == 'lambdef':
                stmt = name
            return context.py__getattribute__(
                name,
                position=stmt.start_pos,
                search_global=True, is_goto=True
            )

    def create_context(self, base_context, node, node_is_context=False, node_is_object=False):
        def parent_scope(node):
            while True:
                node = node.parent

                if parser_utils.is_scope(node):
                    return node
                elif node.type in ('argument', 'testlist_comp'):
                    if node.children[1].type == 'comp_for':
                        return node.children[1]
                elif node.type == 'dictorsetmaker':
                    for n in node.children[1:4]:
                        # In dictionaries it can be pretty much anything.
                        if n.type == 'comp_for':
                            return n

        def from_scope_node(scope_node, child_is_funcdef=None, is_nested=True, node_is_object=False):
            if scope_node == base_node:
                return base_context

            is_funcdef = scope_node.type in ('funcdef', 'lambdef')
            parent_scope = parser_utils.get_parent_scope(scope_node)
            parent_context = from_scope_node(parent_scope, child_is_funcdef=is_funcdef)

            if is_funcdef:
                if isinstance(parent_context, AnonymousInstance):
                    func = BoundMethod(
                        self, parent_context, parent_context.class_context,
                        parent_context.parent_context, scope_node
                    )
                else:
                    func = FunctionContext(
                        self,
                        parent_context,
                        scope_node
                    )
                if is_nested and not node_is_object:
                    return func.get_function_execution()
                return func
            elif scope_node.type == 'classdef':
                class_context = ClassContext(self, parent_context, scope_node)
                if child_is_funcdef:
                    # anonymous instance
                    return AnonymousInstance(self, parent_context, class_context)
                else:
                    return class_context
            elif scope_node.type == 'comp_for':
                if node.start_pos >= scope_node.children[-1].start_pos:
                    return parent_context
                return CompForContext.from_comp_for(parent_context, scope_node)
            raise Exception("There's a scope that was not managed.")

        base_node = base_context.tree_node

        if node_is_context and parser_utils.is_scope(node):
            scope_node = node
        else:
            if node.parent.type in ('funcdef', 'classdef') and node.parent.name == node:
                # When we're on class/function names/leafs that define the
                # object itself and not its contents.
                node = node.parent
            scope_node = parent_scope(node)
        return from_scope_node(scope_node, is_nested=True, node_is_object=node_is_object)

    def parse_and_get_code(self, code=None, path=None, **kwargs):
        if self.allow_different_encoding:
            if code is None:
                with open(path, 'rb') as f:
                    code = f.read()
            code = python_bytes_to_unicode(code, errors='replace')

        return self.grammar.parse(code=code, path=path, **kwargs), code

    def parse(self, *args, **kwargs):
        return self.parse_and_get_code(*args, **kwargs)[0]
