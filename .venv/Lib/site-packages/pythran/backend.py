'''
This module contains all pythran backends.
    * Cxx dumps the AST into C++ code
    * Python dumps the AST into Python code
'''

from pythran.analyses import LocalNodeDeclarations, GlobalDeclarations, Scope
from pythran.analyses import YieldPoints, IsAssigned, ASTMatcher, AST_any
from pythran.analyses import RangeValues, PureExpressions, Dependencies
from pythran.analyses import Immediates, Ancestors, StrictAliases
from pythran.config import cfg
from pythran.cxxgen import Template, Include, Namespace, CompilationUnit
from pythran.cxxgen import Statement, Block, AnnotatedStatement, Typedef, Label
from pythran.cxxgen import Value, FunctionDeclaration, EmptyStatement, Nop
from pythran.cxxgen import FunctionBody, Line, ReturnStatement, Struct, Assign
from pythran.cxxgen import For, While, TryExcept, ExceptHandler, If, AutoFor
from pythran.cxxgen import StatementWithComments, InstrumentedStatement
from pythran.openmp import OMPDirective
from pythran.passmanager import Backend
from pythran.syntax import PythranSyntaxError
from pythran.tables import operator_to_lambda, update_operator_to_lambda
from pythran.tables import pythran_ward, attributes as attributes_table
from pythran.types.conversion import PYTYPE_TO_CTYPE_TABLE, TYPE_TO_SUFFIX
from pythran.types.types import Types
from pythran.utils import attr_to_path, pushpop, cxxid, isstr, isnum
from pythran.utils import isextslice, ispowi, quote_cxxstring
from pythran import metadata, unparse

from math import isnan, isinf
import gast as ast
from functools import reduce
import io

class Python(Backend):
    '''
    Produces a Python representation of the AST.

    >>> import gast as ast, pythran.passmanager as passmanager
    >>> node = ast.parse("print('hello world')")
    >>> pm = passmanager.PassManager('test')
    >>> print(pm.dump(Python, node))
    print('hello world')
    '''

    ResultType = str

    def visit(self, node):
        output = io.StringIO()
        unparse.Unparser(node, output)
        self.result = output.getvalue()


def templatize(node, types, default_types=None):
    if not default_types:
        default_types = [None] * len(types)
    if types:
        return Template(
            ["typename {0} {1}".format(t, "= {0}".format(d) if d else "")
             for t, d in zip(types, default_types)],
            node)
    else:
        return node


def cxx_loop(visit):
    """
    Decorator for loop node (For and While) to handle "else" branching.

    Decorated node will save flags for a goto statement used instead of usual
    break and add this flag at the end of the else statements.

    Examples
    --------
    >> for i in range(12):
    >>     if i == 5:
    >>         break
    >> else:
    >>     ... some code ...

    Becomes

    >> for(type i : range(12))
    >>     if(i==5)
    >>         goto __no_breaking0;
    >> ... some code ...
    >> __no_breaking0;
    """
    def loop_visitor(self, node):
        """
        New decorate function.

        It push the breaking flag, run the visitor and add "else" statements.
        """
        if not node.orelse:
            with pushpop(self.break_handlers, None):
                res = visit(self, node)
            return res

        break_handler = "__no_breaking{0}".format(len(self.break_handlers))
        with pushpop(self.break_handlers, break_handler):
            res = visit(self, node)

        # handle the body of the for loop
        orelse = [self.visit(stmt) for stmt in node.orelse]
        if break_handler in self.used_break:
            orelse_label = [Label(break_handler)]
        else:
            orelse_label = []
        skip = [node.target.id] if isinstance(node, ast.For) else []
        return self.process_locals(node, Block([res] + orelse + orelse_label),
                                   *skip)
    return loop_visitor


class CachedTypeVisitor:

    def __init__(self, other=None):
        if other is None:
            self.mapping = dict()
            self.typeid = dict()
            self.combined = dict()
        else:
            self.mapping = other.mapping.copy()
            self.typeid = other.typeid.copy()
            self.combined = other.combined.copy()

    def __call__(self, node):
        if node not in self.mapping:
            t = node.generate(self)
            if node not in self.mapping:
                # Always re-evaluate LType as their evaluation depends on the
                # callers (due to the recursion clause)
                if type(node).__name__ == 'LType':
                    return t
                else:
                    if t not in self.typeid:
                        self.typeid[t] = len(self.typeid)
                    self.mapping[node] = (t, self.typeid[t])

        return "__type{0}".format(self.mapping[node][1])

    def typedefs(self):
        kv = sorted(set(self.mapping.values()), key=lambda x: x[1])
        L = list()
        for k, v in kv:
            typename = "__type" + str(v)
            L.append(Typedef(Value(k, typename)))
        return L


def make_default(d):
    return "= {0}".format(d) if d else ""


def make_function_declaration(self, node, rtype, name, ftypes, fargs,
                              defaults=None, attributes=None):
    if defaults is None:
        defaults = [None] * len(ftypes)
    if attributes is None:
        attributes = []

    arguments = list()
    first_default = len(node.args.args) - len(node.args.defaults)
    for i, (t, a, d) in enumerate(zip(ftypes, fargs, defaults)):
        argument = Value(t, "{0}{1}".format(a, make_default(d)))
        arguments.append(argument)
    return FunctionDeclaration(Value(rtype, name), arguments, *attributes)


def make_const_function_declaration(self, node, rtype, name, ftypes, fargs,
                                    defaults=None):
    return make_function_declaration(self, node, rtype, name, ftypes, fargs,
                                     defaults, ["const"])


class CxxFunction(ast.NodeVisitor):
    '''
    Attributes
    ----------
    ldecls : {str}
        set of local declarations.
    break_handler : [str]
        It contains flags for goto statements to jump on break in case of
        orelse statement in loop. None means there are no orelse statement so
        no jump are requiered.
        (else in loop means : don't execute if loop is terminated with a break)
    '''

    def __init__(self, parent):
        """ Basic initialiser gathering analysis informations. """
        self.parent = parent
        self.break_handlers = []
        self.used_break = set()
        self.ldecls = None
        self.openmp_deps = set()
        self.unique_counter = 0
        if not (cfg.getboolean('backend', 'annotate') and
                self.passmanager.code):
            self.add_line_info = self.skip_line_info
        else:
            self.lines = self.passmanager.code.split('\n')

    def __getattr__(self, attr):
        return getattr(self.parent, attr)

    def unique(self):
        self.unique_counter += 1
        return self.unique_counter

    # local declaration processing
    def process_locals(self, node, node_visited, *skipped):
        """
        Declare variable local to node and insert declaration before.

        Not possible for function yielding values.
        """
        local_vars = self.scope[node].difference(skipped)
        local_vars = local_vars.difference(self.openmp_deps)
        if not local_vars:
            return node_visited  # no processing

        locals_visited = []
        for varname in sorted(local_vars):
            vartype = self.typeof(varname)
            decl = Statement("{} {}".format(vartype, varname))
            locals_visited.append(decl)
        self.ldecls.difference_update(local_vars)
        return Block(locals_visited + [node_visited])

    def visit_OMPDirective(self, node):
        self.openmp_deps.update(d.id for d in node.private_deps)
        self.openmp_deps.update(d.id for d in node.shared_deps)

    def add_line_info(self, node, cxx_node):
        if not isinstance(node, ast.stmt):
            return cxx_node
        line = self.lines[node.lineno - 1].rstrip()
        if isinstance(node, ast.FunctionDef):
            head, tail = cxx_node
            return head, [StatementWithComments(t, line) for t in tail]
        if cfg.get('backend', 'annotation_kind') == 'lineno':
            return InstrumentedStatement(cxx_node,
                    'pythran_trace_lineno({});'.format(node.lineno))
        else:
            return StatementWithComments(cxx_node, line)

    def skip_line_info(self, node, cxx_node):
        return cxx_node

    def visit(self, node):
        metadata.visit(self, node)
        result = super(CxxFunction, self).visit(node)
        return self.add_line_info(node, result)

    def process_omp_attachements(self, node, stmt, index=None):
        """
        Add OpenMP pragma on the correct stmt in the correct order.

        stmt may be a list. On this case, index have to be specify to add
        OpenMP on the correct statement.
        """
        omp_directives = metadata.get(node, OMPDirective)
        if omp_directives:
            directives = list()
            for directive in omp_directives:
                directive.deps = [self.visit(dep) for dep in directive.deps]
                directives.append(directive)
            if index is None:
                stmt = AnnotatedStatement(stmt, directives)
            else:
                stmt[index] = AnnotatedStatement(stmt[index], directives)
        return stmt

    def typeof(self, node):
        if isinstance(node, str):
            return self.typeof(self.local_names[node])
        elif isinstance(node, ast.AST):
            return self.lctx(self.types[node])
        else:
            return self.lctx(node)

    def prepare_functiondef_context(self, node):
        # prepare context and visit function body
        fargs = node.args.args

        formal_args = [cxxid(arg.id) for arg in fargs]
        formal_types = ["argument_type" + str(i) for i in range(len(fargs))]

        local_decls = set(self.gather(LocalNodeDeclarations, node))
        self.local_names = {sym.id: sym for sym in local_decls}
        self.local_names.update({arg.id: arg for arg in fargs})

        self.lctx = CachedTypeVisitor()

        self.ldecls = {n.id for n in local_decls}
        body = [self.visit(stmt) for stmt in node.body]
        return body, formal_types, formal_args

    def prepare_types(self, node):
        # compute arg dump
        dflt_v = [self.visit(n) for n in node.args.defaults]
        dflt_argv = (
            [None] * (len(node.args.args) - len(node.args.defaults)) +
            dflt_v)
        dflt_argt = (
            [None] * (len(node.args.args) - len(node.args.defaults)) +
            ["decltype({})".format(v) for v in dflt_v])

        # compute type dump
        result_type = self.types[node][0]

        callable_type = Typedef(Value("void", "callable"))
        pure_type = (Typedef(Value("void", "pure"))
                     if node in self.pure_expressions else EmptyStatement())

        return dflt_argv, dflt_argt, result_type, callable_type, pure_type

    # stmt
    def visit_FunctionDef(self, node):

        self.fname = cxxid(node.name)
        tmp = self.prepare_functiondef_context(node)
        operator_body, formal_types, formal_args = tmp

        tmp = self.prepare_types(node)
        dflt_argv, dflt_argt, result_type, callable_type, pure_type = tmp

        # a function has a call operator to be called
        # and a default constructor to create instances
        fscope = "type{0}::".format("<{0}>".format(", ".join(formal_types))
                                    if formal_types
                                    else "")
        ffscope = "{0}::{1}".format(self.fname, fscope)

        operator_declaration = [
            templatize(
                make_const_function_declaration(
                    self, node,
                    "typename {0}result_type".format(fscope),
                    "operator()",
                    formal_types, formal_args, dflt_argv),
                formal_types,
                dflt_argt),
            EmptyStatement()
            ]
        operator_signature = make_const_function_declaration(
            self, node,
            "typename {0}result_type".format(ffscope),
            "{0}::operator()".format(self.fname),
            formal_types, formal_args)

        operator_local_declarations = (
            [Statement("{0} {1}".format(
                self.lctx(self.types[self.local_names[k]]), cxxid(k)))
             for k in sorted(self.ldecls)]
        )
        dependent_typedefs = self.lctx.typedefs()
        operator_definition = FunctionBody(
            templatize(operator_signature, formal_types),
            Block(dependent_typedefs +
                  operator_local_declarations +
                  operator_body)
            )

        ctx = CachedTypeVisitor()
        extra_typedefs = (
            [Typedef(Value(ctx(t), t.name))
             for t in self.types[node][1]] +
            [Typedef(Value(
                ctx(result_type),
                "result_type"))]
        )
        extra_typedefs = ctx.typedefs() + extra_typedefs
        return_declaration = [
            templatize(
                Struct("type", extra_typedefs),
                formal_types,
                dflt_argt
                )
            ]
        topstruct = Struct(self.fname,
                           [callable_type, pure_type] +
                           return_declaration +
                           operator_declaration)

        topstruct = self.process_omp_attachements(node, topstruct)
        return [topstruct], [operator_definition]

    def visit_Return(self, node):
        value = self.visit(node.value)
        if metadata.get(node, metadata.StaticReturn):
            # don't rely on auto because we want to make sure there's no
            # conversion each time we return
            # this happens for variant because the variant param
            # order may differ from the init order (because of the way we
            # do type inference
            rtype = "typename {}::type::result_type".format(self.fname)
            stmt = Block([Assign("static %s tmp_global" % rtype, value),
                          ReturnStatement("tmp_global")])
        else:
            stmt = ReturnStatement(value)
        return self.process_omp_attachements(node, stmt)

    def visit_Delete(self, _):
        return Nop()  # nothing to do in there

    def visit_Assign(self, node):
        """
        Create Assign node for final Cxx representation.

        It tries to handle multi assignment like:

        >> a = b = c = 2

        If only one local variable is assigned, typing is added:

        >> int a = 2;

        TODO: Handle case of multi-assignement for some local variables.

        Finally, process OpenMP clause like #pragma omp atomic
        """
        targets = node.targets if isinstance(node, ast.Assign) else (node.target,)
        if not all(isinstance(n, (ast.Name, ast.Subscript))
                   for n in targets):
            raise PythranSyntaxError(
                "Must assign to an identifier or a subscript",
                node)
        if not node.value:
            return self.visit_Pass(node)
        value = self.visit(node.value)
        stargets = [self.visit(t) for t in targets]
        alltargets = "= ".join(stargets)
        islocal = (len(stargets) == 1 and
                   isinstance(targets[0], ast.Name) and
                   targets[0].id in self.scope[node] and
                   targets[0].id not in self.openmp_deps)
        if islocal:
            # remove this decls from local decls
            self.ldecls.difference_update(t.id for t in targets)
            # add a local declaration
            if self.types[targets[0]].iscombined():
                alltargets = '{} {}'.format(self.typeof(targets[0]),
                                            alltargets)
            elif isinstance(self.types[targets[0]],
                            self.types.builder.Assignable):
                alltargets = '{} {}'.format(
                    self.types.builder.AssignableNoEscape(
                        self.types.builder.NamedType(
                            'decltype({})'.format(value))).sgenerate(),
                    alltargets)
            else:
                assert isinstance(self.types[targets[0]],
                                  (self.types.builder.Lazy,
                                   self.types.builder.NamedType,
                                   self.types.builder.ListType))
                alltargets = '{} {}'.format(
                    self.types.builder.Lazy(
                        self.types.builder.NamedType(
                            'decltype({})'.format(value))).sgenerate(),
                    alltargets)
        stmt = Assign(alltargets, value)
        return self.process_omp_attachements(node, stmt)

    visit_AnnAssign = visit_Assign

    def visit_AugAssign(self, node):
        value = self.visit(node.value)
        target = self.visit(node.target)
        op = update_operator_to_lambda[type(node.op)]
        stmt = Statement(op(target, value)[1:-1])  # strip spurious parenthesis
        return self.process_omp_attachements(node, stmt)

    def visit_Print(self, node):
        values = [self.visit(n) for n in node.values]
        stmt = Statement("pythonic::builtins::print{0}({1})".format(
            "" if node.nl else "_nonl",
            ", ".join(values))
            )
        return self.process_omp_attachements(node, stmt)

    def is_in_collapse(self, loop, node):
        for ancestor in reversed(self.ancestors[loop]):
            if not isinstance(ancestor, ast.For):
                return False
            for directive in metadata.get(ancestor, OMPDirective):
                if 'collapse' in directive.s:
                    # FIXME: check loop depth and range canonicalization
                    if node not in self.pure_expressions:
                        raise PythranSyntaxError(
                            "not pure expression used as loop target inside a "
                            "collapse clause",
                            loop)
                    return True
        assert False, "unreachable state"

    def gen_for(self, node, target, local_iter, local_iter_decl, loop_body):
        """
        Create For representation on iterator for Cxx generation.

        Examples
        --------
        >> "omp parallel for"
        >> for i in range(10):
        >>     ... do things ...

        Becomes

        >> "omp parallel for shared(__iterX)"
        >> for(decltype(__iterX)::iterator __targetX = __iterX.begin();
               __targetX < __iterX.end(); ++__targetX)
        >>         auto&& i = *__targetX;
        >>     ... do things ...

        It the case of not local variable, typing for `i` disappear and typing
        is removed for iterator in case of yields statement in function.
        """
        # Choose target variable for iterator (which is iterator type)
        local_target = "__target{0}".format(self.unique())
        local_target_decl = self.typeof(
            self.types.builder.IteratorOfType(local_iter_decl))

        islocal = (node.target.id not in self.openmp_deps and
                   node.target.id in self.scope[node] and
                   not hasattr(self, 'yields'))
        # If variable is local to the for body it's a ref to the iterator value
        # type
        if islocal:
            local_type = "auto&&"
            self.ldecls.remove(node.target.id)
        else:
            local_type = ""

        # Assign iterable value
        loop_body_prelude = Statement("{} {}= *{}".format(local_type,
                                                          target,
                                                          local_target))

        # Create the loop
        assign = self.make_assign(local_target_decl, local_target, local_iter)
        loop = For("{}.begin()".format(assign),
                   "{0} < {1}.end()".format(local_target, local_iter),
                   "++{0}".format(local_target),
                   Block([loop_body_prelude, loop_body]))
        return [self.process_omp_attachements(node, loop)]

    def handle_real_loop_comparison(self, args, target, upper_bound):
        """
        Handle comparison for real loops.

        Add the correct comparison operator if possible.
        """
        # order is 1 for increasing loop, -1 for decreasing loop and 0 if it is
        # not known at compile time
        if len(args) <= 2:
            order = 1
        elif isnum(args[2]):
            order = -1 + 2 * (int(args[2].value) > 0)
        elif isnum(args[1]) and isnum(args[0]):
            order = -1 + 2 * (int(args[1].value) > int(args[0].value))
        else:
            order = 0

        comparison = "{} < {}" if order == 1 else "{} > {}"
        comparison = comparison.format(target, upper_bound)
        return comparison

    def gen_c_for(self, node, local_iter, loop_body):
        """
        Create C For representation for Cxx generation.

        Examples
        --------
        >> for i in range(10):
        >>     ... do things ...

        Becomes

        >> for(long i = 0, __targetX = 10; i < __targetX; i += 1)
        >>     ... do things ...

        Or

        >> for i in range(10, 0, -1):
        >>     ... do things ...

        Becomes

        >> for(long i = 10, __targetX = 0; i > __targetX; i += -1)
        >>     ... do things ...


        It the case of not local variable, typing for `i` disappear
        """
        args = node.iter.args
        step = "1L" if len(args) <= 2 else self.visit(args[2])
        if len(args) == 1:
            lower_bound = "0L"
            upper_arg = 0
        else:
            lower_bound = self.visit(args[0])
            upper_arg = 1

        upper_type = iter_type = "long "
        upper_value = self.visit(args[upper_arg])
        if self.is_in_collapse(node, args[upper_arg]):
            upper_bound = upper_value  # compatible with collapse
        else:
            upper_bound = "__target{0}".format(self.unique())

        islocal = (node.target.id not in self.openmp_deps and
                   node.target.id in self.scope[node] and
                   not hasattr(self, 'yields'))
        # If variable is local to the for body keep it local...
        if islocal:
            loop = list()
            self.ldecls.remove(node.target.id)
        else:
            # For yield function, upper_bound is globals.
            iter_type = ""

            # Back one step to keep Python behavior (unless the loop index does
            # not escape)
            if node.target.id in self.scope[node]:
                loop = []
            else:
                loop = [If("{} == {}".format(local_iter, upper_bound),
                        Statement("{} -= {}".format(local_iter, step)))]

        comparison = self.handle_real_loop_comparison(args, local_iter,
                                                      upper_bound)

        forloop = For("{0} {1}={2}".format(iter_type, local_iter, lower_bound),
                      comparison,
                      "{0} += {1}".format(local_iter, step),
                      loop_body)

        loop.insert(0, self.process_omp_attachements(node, forloop))

        # Store upper bound value if needed
        if upper_bound is upper_value:
            header = []
        else:
            assgnt = self.make_assign(upper_type, upper_bound, upper_value)
            header = [Statement(assgnt)]
        return header, loop

    def handle_omp_for(self, node, local_iter):
        """
        Fix OpenMP directives on For loops.

        Add the target as private variable as a new variable may have been
        introduce to handle cxx iterator.

        Also, add the iterator as shared variable as all 'parallel for chunck'
        have to use the same iterator.
        """
        for directive in metadata.get(node, OMPDirective):
            if any(key in directive.s for key in (' parallel ', ' task ')):
                # Eventually add local_iter in a shared clause as iterable is
                # shared in the for loop (for every clause with datasharing)
                directive.s += ' shared({})'
                directive.deps.append(ast.Name(local_iter, ast.Load(),
                                               None, None))
                directive.shared_deps.append(directive.deps[-1])

            target = node.target
            assert isinstance(target, ast.Name)
            hasfor = 'for' in directive.s
            nodefault = 'default' not in directive.s
            noindexref = all(isinstance(x, ast.Name) and
                             x.id != target.id for x in directive.deps)
            if (hasfor and nodefault and noindexref and
                    target.id not in self.scope[node]):
                # Target is private by default in omp but iterator use may
                # introduce an extra variable
                directive.s += ' private({})'
                directive.deps.append(ast.Name(target.id, ast.Load(),
                                               None, None))
                directive.private_deps.append(directive.deps[-1])

    def can_use_autofor(self, node):
        """
        Check if given for Node can use autoFor syntax.

        To use auto_for:
            - iterator should have local scope
            - yield should not be use
            - OpenMP pragma should not be use

        TODO : Yield should block only if it is use in the for loop, not in the
               whole function.
        """
        auto_for = (isinstance(node.target, ast.Name) and
                    node.target.id in self.scope[node] and
                    node.target.id not in self.openmp_deps)
        auto_for &= not metadata.get(node, OMPDirective)
        return auto_for

    def can_use_c_for(self, node):
        """
        Check if a for loop can use classic C syntax.

        To use C syntax:
            - target should not be assign in the loop
            - range should be use as iterator
            - order have to be known at compile time
        """
        assert isinstance(node.target, ast.Name)
        pattern_range = ast.Call(func=ast.Attribute(
            value=ast.Name('builtins', ast.Load(), None, None),
            attr='range', ctx=ast.Load()),
            args=AST_any(), keywords=[])
        is_assigned = set()
        for stmt in node.body:
            is_assigned.update({n.id for n in self.gather(IsAssigned, stmt)})

        match = ASTMatcher(pattern_range).match(node.iter)
        if not match:
            return False
        if node.target.id in is_assigned:
            return False

        args = node.iter.args
        if len(args) < 3:
            return True
        if isnum(args[2]):
            return True
        return False

    def make_assign(self, local_iter_decl, local_iter, iterable):
        return "{0} {1} = {2}".format(local_iter_decl, local_iter, iterable)

    def is_user_function(self, func):
        aliases = self.strict_aliases[func]
        if not aliases:
            return False
        for alias in aliases:
            if not isinstance(alias, ast.FunctionDef):
                return False
            if self.gather(YieldPoints, alias):
                return False
        return True

    @cxx_loop
    def visit_For(self, node):
        """
        Create For representation for Cxx generation.

        Examples
        --------
        >> for i in range(10):
        >>     ... work ...

        Becomes

        >> typename returnable<decltype(builtins.range(10))>::type __iterX
           = builtins.range(10);
        >> ... possible container size reservation ...
        >> for (auto&& i: __iterX)
        >>     ... the work ...

        This function also handle assignment for local variables.

        We can notice that three kind of loop are possible:
        - Normal for loop on iterator
        - Autofor loop.
        - Normal for loop using integer variable iteration
        Kind of loop used depend on OpenMP, yield use and variable scope.
        """
        if not isinstance(node.target, ast.Name):
            raise PythranSyntaxError(
                "Using something other than an identifier as loop target",
                node.target)
        target = self.visit(node.target)

        # Handle the body of the for loop
        loop_body = Block([self.visit(stmt) for stmt in node.body])

        # Declare local variables at the top of the loop body
        if not node.orelse:
            loop_body = self.process_locals(node, loop_body, node.target.id)
        iterable = self.visit(node.iter)

        if self.can_use_c_for(node):
            header, loop = self.gen_c_for(node, target, loop_body)
        else:

            if self.can_use_autofor(node):
                header = []
                self.ldecls.remove(node.target.id)
                autofor = AutoFor(target, iterable, loop_body)
                loop = [self.process_omp_attachements(node, autofor)]
            else:
                # Iterator declaration
                local_iter = "__iter{0}".format(self.unique())
                local_iter_decl = self.types.builder.Assignable(
                    self.types[node.iter])

                self.handle_omp_for(node, local_iter)

                # Assign iterable
                # For C loop, it avoids issues
                # if the upper bound is assigned in the loop
                asgnt = self.make_assign(self.typeof(local_iter_decl),
                                         local_iter, iterable)
                header = [Statement(asgnt)]
                loop = self.gen_for(node, target, local_iter, local_iter_decl,
                                    loop_body)

        # For xxxComprehension, it is replaced by a for loop. In this case,
        # pre-allocate size of container.
        for comp in metadata.get(node, metadata.Comprehension):
            header.append(Statement("pythonic::utils::reserve({0},{1})".format(
                comp.target,
                iterable)))

        return Block(header + loop)

    @cxx_loop
    def visit_While(self, node):
        """
        Create While node for Cxx generation.

        It is a cxx_loop to handle else clause.
        """
        test = self.visit(node.test)
        body = [self.visit(n) for n in node.body]
        stmt = While(test, Block(body))
        return self.process_omp_attachements(node, stmt)

    def visit_Try(self, node):
        body = [self.visit(n) for n in node.body]
        except_ = list()
        for n in node.handlers:
            except_.extend(self.visit(n))
        return TryExcept(Block(body), except_)

    def visit_ExceptHandler(self, node):
        name = self.visit(node.name) if node.name else None
        body = [self.visit(m) for m in node.body]
        if isinstance(node.type, ast.Tuple):
            return [ExceptHandler(p.attr, Block(body), name)
                    for p in node.type.elts]
        else:
            return [ExceptHandler(
                node.type and node.type.attr,
                Block(body),
                name)]

    def visit_If(self, node):
        test = self.visit(node.test)
        body = [self.visit(n) for n in node.body]
        orelse = [self.visit(n) for n in node.orelse]
        # compound statement required for some OpenMP Directives
        if isnum(node.test) and node.test.value == 1:
            stmt = Block(body)
        else:
            stmt = If(test, Block(body), Block(orelse) if orelse else None)
        return self.process_locals(node,
                                   self.process_omp_attachements(node, stmt))

    def visit_Raise(self, node):
        exc = node.exc and self.visit(node.exc)
        return Statement("throw {0}".format(exc or ""))

    def visit_Assert(self, node):
        params = [self.visit(node.test), node.msg and self.visit(node.msg)]
        sparams = ", ".join(_f for _f in params if _f)
        return Statement("pythonic::pythran_assert({0})".format(sparams))

    def visit_Import(self, _):
        return Nop()  # everything is already #included

    def visit_ImportFrom(self, _):
        assert False, "should be filtered out by the expand_import pass"

    def visit_Expr(self, node):
        stmt = Statement(self.visit(node.value))
        return self.process_locals(node,
                                   self.process_omp_attachements(node, stmt))

    def visit_Pass(self, node):
        stmt = EmptyStatement()
        return self.process_omp_attachements(node, stmt)

    def visit_Break(self, _):
        """
        Generate break statement in most case and goto for orelse clause.

        See Also : cxx_loop
        """
        if self.break_handlers and self.break_handlers[-1]:
            self.used_break.add(self.break_handlers[-1])
            return Statement("goto {0}".format(self.break_handlers[-1]))
        else:
            return Statement("break")

    def visit_Continue(self, _):
        return Statement("continue")

    # expr
    def visit_BoolOp(self, node):
        values = [self.visit(value) for value in node.values]
        op = operator_to_lambda[type(node.op)]
        return reduce(op, values)

    def visit_BinOp(self, node):
        left = self.visit(node.left)
        right = self.visit(node.right)
        # special case pow for positive integral exponent
        if ispowi(node):
            right = 'std::integral_constant<long, {}>{{}}'.format(
                node.right.value)
        if isstr(node.left):
            left = "pythonic::types::str({})".format(left)
        elif isstr(node.right):
            right = "pythonic::types::str({})".format(right)
        return operator_to_lambda[type(node.op)](left, right)

    def visit_UnaryOp(self, node):
        operand = self.visit(node.operand)
        return operator_to_lambda[type(node.op)](operand)

    def visit_IfExp(self, node):
        test = self.visit(node.test)
        body = self.visit(node.body)
        orelse = self.visit(node.orelse)
        return (
            "(pythonic::builtins::functor::bool_{{}}({0}) "
            "? typename __combined<decltype({1}), decltype({2})>::type({1}) "
            ": typename __combined<decltype({1}), decltype({2})>::type({2}))"
        ).format(test, body, orelse)

    def visit_List(self, node):
        if not node.elts:  # empty list
            return '{}(pythonic::types::empty_list())'.format(self.typeof(node))
        else:
            node_type = self.types.builder.Assignable(self.types[node])
            elts = [self.visit(n) for n in node.elts]
            return "{0}({{{1}}})".format(
                self.typeof(node_type),
                ", ".join("static_cast<typename {}::value_type>({})"
                          .format(self.typeof(node_type), elt) for elt in elts))

    def visit_Set(self, node):
        if not node.elts:  # empty set
            return '{}(pythonic::types::empty_set())'.format(self.typeof(node))
        else:
            elts = [self.visit(n) for n in node.elts]
            node_type = self.types.builder.Assignable(self.types[node])
            return "{0}({{{1}}})".format(
                self.typeof(node_type),
                ", ".join("static_cast<typename {}::value_type>({})"
                          .format(self.typeof(node_type), elt) for elt in elts))

    def visit_Dict(self, node):
        if not node.keys:  # empty dict
            return '{}(pythonic::types::empty_dict())'.format(self.typeof(node))
        else:
            keys = [self.visit(n) for n in node.keys]
            values = [self.visit(n) for n in node.values]
            return "{0}{{{{{1}}}}}".format(
                self.typeof(self.types.builder.Assignable(self.types[node])),
                ", ".join("{{ {0}, {1} }}".format(k, v)
                          for k, v in zip(keys, values)))

    def visit_Tuple(self, node):
        elts = [self.visit(elt) for elt in node.elts]
        tuple_type = self.types[node]
        result = "pythonic::types::make_tuple({0})".format(", ".join(elts))
        if isinstance(tuple_type, self.types.builder.CombinedTypes):
            return '({}){}'.format(self.typeof(tuple_type), result)
        else:
            return result

    def visit_Compare(self, node):
        left = self.visit(node.left)
        ops = [operator_to_lambda[type(n)] for n in node.ops]
        comparators = [self.visit(n) for n in node.comparators]
        all_cmps = zip([left] + comparators[:-1], ops, comparators)
        return " and ".join(op(x, y) for x, op, y in all_cmps)

    def visit_Call(self, node):
        args = [self.visit(n) for n in node.args]
        func = self.visit(node.func)
        # special hook for getattr, as we cannot represent it in C++
        if func == 'pythonic::builtins::functor::getattr{}':
            attrname = node.args[1].value
            fmt = 'pythonic::builtins::getattr({}{{}}, {})'
            attr = 'pythonic::types::attr::' + attrname.upper()
            if attributes_table[attrname][1].isstatic() and node in self.immediates:
                # ugly hack to ensure constexprness of the call
                arg = '(decltype(&{}))nullptr'.format(args[0])
            else:
                arg = args[0]
            result = fmt.format(attr, arg)
        # Avoid passing scalars by ref as it prevents some C++ optimization.
        # pythonic::types::call (tries to) handle that gracefully.
        elif args and self.is_user_function(node.func):
            result = "pythonic::types::call({})".format(", ".join([func] + args))
        else:
            result = "{}({})".format(func, ", ".join(args))

        # When we have extra type information to inject as a cast
        if isinstance(self.types.get(node), self.types.builder.CombinedTypes):
            return '({}){}'.format(self.typeof(node), result)
        else:
            return result

    def visit_Constant(self, node):
        if node.value is None:
            ret = 'pythonic::builtins::None'
        elif isinstance(node.value, bool):
            ret = str(node.value).lower()
        elif isinstance(node.value, bytes):
            quoted = "".join('\\' + hex(b)[1:] for b in node.value)
            # FIXME: using str type as backend
            if len(node.value) == 1:
                quoted = quoted.replace("'", r"\'")
                ret = 'pythonic::types::chr(\'' + quoted + '\')'
            else:
                ret = 'pythonic::types::str("' + quoted + '")'
        elif isinstance(node.value, str):
            quoted = quote_cxxstring(node.value)
            if len(node.value) == 1:
                quoted = quoted.replace("'", r"\'")
                ret = 'pythonic::types::chr(\'' + quoted + '\')'
            else:
                ret = 'pythonic::types::str("' + quoted + '")'
        elif isinstance(node.value, complex):
            ret = "{0}({1}, {2})".format(
                PYTYPE_TO_CTYPE_TABLE[complex],
                node.value.real,
                node.value.imag)
        elif isnan(node.value):
            ret = 'pythonic::numpy::nan'
        elif isinf(node.value):
            ret = ('+' if node.value >= 0 else '-') + 'pythonic::numpy::inf'
        else:
            ret = repr(node.value) + TYPE_TO_SUFFIX.get(type(node.value), "")
        if node in self.immediates:
            if isinstance(node.value, int):
                return "std::integral_constant<%s, %s>{}" % (
                    PYTYPE_TO_CTYPE_TABLE[type(node.value)], str(node.value).lower())
            if isinstance(node.value, str):
                assert len(node.value) == 1
                return "std::integral_constant<char, '%s'>{}" % node.value
            raise PythranSyntaxError("Unsupported immediate type", node)
        return ret

    def visit_Attribute(self, node):
        obj, path = attr_to_path(node)
        sattr = '::'.join(map(cxxid, path))
        if not obj.isliteral():
            sattr += '{}'
        return sattr

    def all_positive(self, node):
        if isinstance(node, ast.Tuple):
            return all(self.range_values[elt].low >= 0
                       for elt in node.elts)
        return self.range_values[node].low >= 0

    def stores_to(self, node):
        ancestors = self.ancestors[node] + (node,)
        stmt_indices = [i for i, n in enumerate(ancestors)
                        if isinstance(n, (ast.Assign, ast.For))]
        if not stmt_indices:
            return True

        stmt_index = stmt_indices[-1]

        if isinstance(ancestors[stmt_index], ast.Assign):
            return ancestors[stmt_index + 1] is ancestors[stmt_index].value
        else:
            return ancestors[stmt_index + 1] is not ancestors[stmt_index].target

    def visit_Subscript(self, node):
        value = self.visit(node.value)
        if self.stores_to(node):
            value = 'pythonic::types::as_const({})'.format(value)

        # we cannot overload the [] operator in that case
        if isstr(node.value):
            value = 'pythonic::types::str({})'.format(value)
        # positive static index case
        if (isnum(node.slice) and
                (node.slice.value >= 0) and
                isinstance(node.slice.value, int)):
            return "std::get<{0}>({1})".format(node.slice.value, value)
        # positive indexing case
        elif self.all_positive(node.slice):
            slice_ = self.visit(node.slice)
            return "{1}.fast({0})".format(slice_, value)
        # extended slice case
        elif isextslice(node.slice):
            slices = [self.visit(elt) for elt in node.slice.elts]
            return "{1}({0})".format(','.join(slices), value)
        # standard case
        else:
            slice_ = self.visit(node.slice)
            return "{1}[{0}]".format(slice_, value)

    def visit_Name(self, node):
        if node.id in self.local_names:
            return cxxid(node.id)
        elif node.id in self.global_declarations:
            return "{0}()".format(cxxid(node.id))
        else:
            return cxxid(node.id)

    # other

    def visit_Slice(self, node):
        args = []
        for field in ('lower', 'upper', 'step'):
            nfield = getattr(node, field)
            arg = (self.visit(nfield) if nfield
                   else 'pythonic::builtins::None')
            args.append(arg)

        nstep = node.step
        if nstep is None or (isnum(nstep) and nstep.value > 0):
            if nstep is None or nstep.value == 1:
                if self.all_positive(node.lower) and self.all_positive(node.upper):
                    builder = "pythonic::types::fast_contiguous_slice({0},{1})"
                else:
                    builder = "pythonic::types::contiguous_slice({0},{1})"
                step = 1
            else:
                builder = "pythonic::types::cstride_slice<{2}>({0},{1})"
                step = nstep.value

            return builder.format(args[0], args[1], step)
        else:
            return "pythonic::types::slice({},{},{})".format(*args)


class CxxGenerator(CxxFunction):

    # recover previous generator state
    StateHolder = "__generator_state"
    StateValue = "__generator_value"
    # flags the last statement of a generator
    FinalStatement = "that_is_all_folks"
    # local declaration processing

    def process_locals(self, node, node_visited, *skipped):
        return node_visited  # no processing

    def prepare_functiondef_context(self, node):
        self.extra_declarations = []

        # 0 is used as initial_state, thus the +1
        self.yields = {k: (1 + v, "yield_point{0}".format(1 + v)) for (v, k) in
                       enumerate(self.gather(YieldPoints, node))}
        return super(CxxGenerator, self).prepare_functiondef_context(node)

    # stmt
    def visit_FunctionDef(self, node):
        self.returns = False
        tmp = self.prepare_functiondef_context(node)
        operator_body, formal_types, formal_args = tmp

        tmp = self.prepare_types(node)
        dflt_argv, dflt_argt, result_type, callable_type, pure_type = tmp

        # a generator has a call operator that returns the iterator
        next_name = "__generator__{0}".format(cxxid(node.name))
        instanciated_next_name = "{0}{1}".format(
            next_name,
            "<{0}>".format(", ".join(formal_types)) if formal_types else "")

        if self.returns:
            operator_body.append(Label(CxxGenerator.FinalStatement))
        operator_body.append(Statement("return result_type()"))

        next_declaration = [
            FunctionDeclaration(Value("result_type", "next"), []),
            EmptyStatement()]  # empty statement to force a comma ...

        # the constructors
        next_constructors = [
            FunctionBody(
                FunctionDeclaration(Value("", next_name), []),
                Line(': pythonic::yielder() {}')
                )]
        if formal_types:
            # If all parameters have a default value, we don't need default
            # constructor
            if dflt_argv and all(dflt_argv):
                next_constructors = list()
            next_constructors.append(FunctionBody(
                make_function_declaration(self, node, "", next_name,
                                          formal_types, formal_args,
                                          dflt_argv),
                Line(": {0} {{ }}".format(
                    ", ".join(["pythonic::yielder()"] +
                              ["{0}({0})".format(arg)
                               for arg in formal_args])))
                ))

        next_iterator = [
            FunctionBody(
                FunctionDeclaration(Value("void", "operator++"), []),
                Block([Statement("next()")])),
            FunctionBody(
                FunctionDeclaration(
                    Value("typename {0}::result_type".format(
                        instanciated_next_name),
                        "operator*"),
                    [], "const"),
                Block([
                    ReturnStatement(
                        CxxGenerator.StateValue)])),
            FunctionBody(
                FunctionDeclaration(
                    Value("pythonic::types::generator_iterator<{0}>"
                          .format(next_name),
                          "begin"),
                    []),
                Block([Statement("next()"),
                       ReturnStatement(
                           "pythonic::types::generator_iterator<{0}>"
                           "(*this)".format(next_name))])),
            FunctionBody(
                FunctionDeclaration(
                    Value("pythonic::types::generator_iterator<{0}>"
                          .format(next_name),
                          "end"),
                    []),
                Block([ReturnStatement(
                    "pythonic::types::generator_iterator<{0}>()"
                    .format(next_name))]))
            ]
        next_signature = templatize(
            FunctionDeclaration(
                Value(
                    "typename {0}::result_type".format(
                        instanciated_next_name),
                    "{0}::next".format(instanciated_next_name)),
                []),
            formal_types)

        next_body = operator_body
        # the dispatch table at the entry point
        next_body.insert(0, Statement("switch({0}) {{ {1} }}".format(
            CxxGenerator.StateHolder,
            " ".join("case {0}: goto {1};".format(num, where)
                     for (num, where) in sorted(
                         self.yields.values(),
                         key=lambda x: x[0])))))

        ctx = CachedTypeVisitor(self.lctx)
        next_members = ([Statement("{0} {1}".format(ft, fa))
                         for (ft, fa) in zip(formal_types, formal_args)] +
                        [Statement("{0} {1}".format(
                            ctx(self.types[self.local_names[k]]),
                            k))
                         for k in sorted(self.ldecls)] +
                        [Statement("{0} {1}".format(v, k))
                         for k, v in self.extra_declarations] +
                        [Statement(
                            "typename {0}::result_type {1}".format(
                                instanciated_next_name,
                                CxxGenerator.StateValue))])

        extern_typedefs = [Typedef(Value(ctx(t), t.name))
                           for t in self.types[node][1]]
        iterator_typedef = [
            Typedef(
                Value("pythonic::types::generator_iterator<{0}>".format(
                    "{0}<{1}>".format(next_name, ", ".join(formal_types))
                    if formal_types else next_name),
                    "iterator")),
            Typedef(Value(ctx(result_type),
                          "value_type"))]
        result_typedef = [
            Typedef(Value(ctx(result_type), "result_type"))]
        extra_typedefs = (ctx.typedefs() +
                          extern_typedefs +
                          iterator_typedef +
                          result_typedef)

        next_struct = templatize(
            Struct(next_name,
                   extra_typedefs +
                   next_members +
                   next_constructors +
                   next_iterator +
                   next_declaration, "pythonic::yielder"),
            formal_types)
        next_definition = FunctionBody(next_signature, Block(next_body))

        operator_declaration = [
            templatize(
                make_const_function_declaration(
                    self, node, instanciated_next_name,
                    "operator()",
                    formal_types, formal_args, dflt_argv),
                formal_types,
                dflt_argt),
            EmptyStatement()]
        operator_signature = make_const_function_declaration(
            self, node, instanciated_next_name,
            "{0}::operator()".format(cxxid(node.name)),
            formal_types, formal_args)
        operator_definition = FunctionBody(
            templatize(operator_signature, formal_types),
            Block([ReturnStatement("{0}({1})".format(
                instanciated_next_name,
                ", ".join(formal_args)))])
            )

        topstruct_type = templatize(
            Struct("type", extra_typedefs),
            formal_types)
        topstruct = Struct(
            cxxid(node.name),
            [topstruct_type, callable_type, pure_type] +
            operator_declaration)

        return [next_struct, topstruct], [next_definition, operator_definition]

    def visit_Return(self, node):
        self.returns = True
        return Block([Statement("{0} = -1".format(CxxGenerator.StateHolder)),
                      Statement("goto {0}".format(CxxGenerator.FinalStatement))
                      ])

    def visit_Yield(self, node):
        num, label = self.yields[node]
        return "".join(n for n in Block([
            Assign(CxxGenerator.StateHolder, num),
            ReturnStatement("{0} = {1}".format(CxxGenerator.StateValue,
                                               self.visit(node.value))),
            Statement("{0}:".format(label))
            ]).generate())

    def visit_Assign(self, node):
        if not node.value:
            return self.visit_Pass(node)
        value = self.visit(node.value)
        targets = node.targets if isinstance(node, ast.Assign) else (node.target,)
        targets = [self.visit(t) for t in targets]
        alltargets = "= ".join(targets)
        stmt = Assign(alltargets, value)
        return self.process_omp_attachements(node, stmt)

    visit_AnnAssign = visit_Assign

    def can_use_autofor(self, node):
        """
        TODO : Yield should block only if it is use in the for loop, not in the
               whole function.
        """
        return False

    def make_assign(self, local_iter_decl, local_iter, iterable):
        # For yield function, iterable is globals.
        self.extra_declarations.append((local_iter, local_iter_decl,))
        return super(CxxGenerator, self).make_assign("", local_iter, iterable)


class Cxx(Backend[Dependencies, GlobalDeclarations, Types, Scope, RangeValues,
                  PureExpressions, Immediates, Ancestors, StrictAliases]):

    """
    Produces a C++ representation of the AST.

    >>> import gast as ast, pythran.passmanager as passmanager, os
    >>> node = ast.parse("def foo(): return 'hello world'")
    >>> pm = passmanager.PassManager('test')
    >>> r = pm.dump(Cxx, node)
    >>> print(str(r).replace(os.sep, '/'))
    #include <pythonic/include/types/str.hpp>
    #include <pythonic/types/str.hpp>
    namespace 
    {
      namespace __pythran_test
      {
        struct foo
        {
          typedef void callable;
          typedef void pure;
          struct type
          {
            typedef pythonic::types::str __type0;
            typedef typename pythonic::returnable<__type0>::type __type1;
            typedef __type1 result_type;
          }  ;
          inline
          typename type::result_type operator()() const;
          ;
        }  ;
        inline
        typename foo::type::result_type foo::operator()() const
        {
          return pythonic::types::str("hello world");
        }
      }
    }
    """
    ResultType = type(None)

    # mod
    def visit_Module(self, node):
        """ Build a compilation unit. """
        if cfg.getboolean('backend', 'annotate'):
            node = ast.fix_missing_locations(node)
        # build all types
        header_deps = sorted(self.dependencies)
        headers = [Include('/'.join(["pythonic", "include"] +
                                    [cxxid(x) for x in t]) + ".hpp")
                   for t in header_deps]
        headers += [Include('/'.join(["pythonic"] + [cxxid(x) for x in t])
                            + ".hpp")
                    for t in header_deps]

        decls_n_defns = list(filter(None, (self.visit(stmt) for stmt in
                                    node.body)))
        decls, defns = zip(*decls_n_defns) if decls_n_defns else ([], [])

        nsbody = [s for ls in decls + defns for s in ls]
        ns = Namespace(pythran_ward + self.passmanager.module_name, nsbody)
        anonymous_ns = Namespace("", [ns])  # force internal linkage
        self.result = CompilationUnit(headers + [anonymous_ns])

    def visit_FunctionDef(self, node):
        yields = self.gather(YieldPoints, node)
        visitor = (CxxGenerator if yields else CxxFunction)(self)
        return visitor.visit(node)
