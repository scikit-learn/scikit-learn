""" Turn compiler.ast structures back into executable python code.

    The unparse method takes a compiler.ast tree and transforms it back into
    valid python code.  It is incomplete and currently only works for
    import statements, function calls, function definitions, assignments, and
    basic expressions.

    Inspired by python-2.5-svn/Demo/parser/unparse.py

    fixme: We may want to move to using _ast trees because the compiler for
           them is about 6 times faster than compiler.compile.
"""
from __future__ import division, absolute_import, print_function

import sys
from compiler.ast import Const, Name, Tuple, Div, Mul, Sub, Add

if sys.version_info[0] >= 3:
    from io import StringIO
else:
    from StringIO import StringIO

def unparse(ast, single_line_functions=False):
    s = StringIO()
    UnparseCompilerAst(ast, s, single_line_functions)
    return s.getvalue().lstrip()

op_precedence = { 'compiler.ast.Power':3, 'compiler.ast.Mul':2, 'compiler.ast.Div':2,
                  'compiler.ast.Add':1, 'compiler.ast.Sub':1 }

class UnparseCompilerAst:
    """ Methods in this class recursively traverse an AST and
        output source code for the abstract syntax; original formatting
        is disregarged.
    """

    #########################################################################
    # object interface.
    #########################################################################

    def __init__(self, tree, file = sys.stdout, single_line_functions=False):
        """ Unparser(tree, file=sys.stdout) -> None.

            Print the source for tree to file.
        """
        self.f = file
        self._single_func = single_line_functions
        self._do_indent = True
        self._indent = 0
        self._dispatch(tree)
        self._write("\n")
        self.f.flush()

    #########################################################################
    # Unparser private interface.
    #########################################################################

    ### format, output, and dispatch methods ################################

    def _fill(self, text = ""):
        "Indent a piece of text, according to the current indentation level"
        if self._do_indent:
            self._write("\n"+"    "*self._indent + text)
        else:
            self._write(text)

    def _write(self, text):
        "Append a piece of text to the current line."
        self.f.write(text)

    def _enter(self):
        "Print ':', and increase the indentation."
        self._write(": ")
        self._indent += 1

    def _leave(self):
        "Decrease the indentation level."
        self._indent -= 1

    def _dispatch(self, tree):
        "_dispatcher function, _dispatching tree type T to method _T."
        if isinstance(tree, list):
            for t in tree:
                self._dispatch(t)
            return
        meth = getattr(self, "_"+tree.__class__.__name__)
        if tree.__class__.__name__ == 'NoneType' and not self._do_indent:
            return
        meth(tree)


    #########################################################################
    # compiler.ast unparsing methods.
    #
    # There should be one method per concrete grammar type. They are
    # organized in alphabetical order.
    #########################################################################

    def _Add(self, t):
        self.__binary_op(t, '+')

    def _And(self, t):
        self._write(" (")
        for i, node in enumerate(t.nodes):
            self._dispatch(node)
            if i != len(t.nodes)-1:
                self._write(") and (")
        self._write(")")

    def _AssAttr(self, t):
        """ Handle assigning an attribute of an object
        """
        self._dispatch(t.expr)
        self._write('.'+t.attrname)

    def _Assign(self, t):
        """ Expression Assignment such as "a = 1".

            This only handles assignment in expressions.  Keyword assignment
            is handled separately.
        """
        self._fill()
        for target in t.nodes:
            self._dispatch(target)
            self._write(" = ")
        self._dispatch(t.expr)
        if not self._do_indent:
            self._write('; ')

    def _AssName(self, t):
        """ Name on left hand side of expression.

            Treat just like a name on the right side of an expression.
        """
        self._Name(t)

    def _AssTuple(self, t):
        """ Tuple on left hand side of an expression.
        """

        # _write each elements, separated by a comma.
        for element in t.nodes[:-1]:
            self._dispatch(element)
            self._write(", ")

        # Handle the last one without writing comma
        last_element = t.nodes[-1]
        self._dispatch(last_element)

    def _AugAssign(self, t):
        """ +=,-=,*=,/=,**=, etc. operations
        """

        self._fill()
        self._dispatch(t.node)
        self._write(' '+t.op+' ')
        self._dispatch(t.expr)
        if not self._do_indent:
            self._write(';')

    def _Bitand(self, t):
        """ Bit and operation.
        """

        for i, node in enumerate(t.nodes):
            self._write("(")
            self._dispatch(node)
            self._write(")")
            if i != len(t.nodes)-1:
                self._write(" & ")

    def _Bitor(self, t):
        """ Bit or operation
        """

        for i, node in enumerate(t.nodes):
            self._write("(")
            self._dispatch(node)
            self._write(")")
            if i != len(t.nodes)-1:
                self._write(" | ")

    def _CallFunc(self, t):
        """ Function call.
        """
        self._dispatch(t.node)
        self._write("(")
        comma = False
        for e in t.args:
            if comma: self._write(", ")
            else: comma = True
            self._dispatch(e)
        if t.star_args:
            if comma: self._write(", ")
            else: comma = True
            self._write("*")
            self._dispatch(t.star_args)
        if t.dstar_args:
            if comma: self._write(", ")
            else: comma = True
            self._write("**")
            self._dispatch(t.dstar_args)
        self._write(")")

    def _Compare(self, t):
        self._dispatch(t.expr)
        for op, expr in t.ops:
            self._write(" " + op + " ")
            self._dispatch(expr)

    def _Const(self, t):
        """ A constant value such as an integer value, 3, or a string, "hello".
        """
        self._dispatch(t.value)

    def _Decorators(self, t):
        """ Handle function decorators (eg. @has_units)
        """
        for node in t.nodes:
            self._dispatch(node)

    def _Dict(self, t):
        self._write("{")
        for  i, (k, v) in enumerate(t.items):
            self._dispatch(k)
            self._write(": ")
            self._dispatch(v)
            if i < len(t.items)-1:
                self._write(", ")
        self._write("}")

    def _Discard(self, t):
        """ Node for when return value is ignored such as in "foo(a)".
        """
        self._fill()
        self._dispatch(t.expr)

    def _Div(self, t):
        self.__binary_op(t, '/')

    def _Ellipsis(self, t):
        self._write("...")

    def _From(self, t):
        """ Handle "from xyz import foo, bar as baz".
        """
        # fixme: Are From and ImportFrom handled differently?
        self._fill("from ")
        self._write(t.modname)
        self._write(" import ")
        for i, (name,asname) in enumerate(t.names):
            if i != 0:
                self._write(", ")
            self._write(name)
            if asname is not None:
                self._write(" as "+asname)

    def _Function(self, t):
        """ Handle function definitions
        """
        if t.decorators is not None:
            self._fill("@")
            self._dispatch(t.decorators)
        self._fill("def "+t.name + "(")
        defaults = [None] * (len(t.argnames) - len(t.defaults)) + list(t.defaults)
        for i, arg in enumerate(zip(t.argnames, defaults)):
            self._write(arg[0])
            if arg[1] is not None:
                self._write('=')
                self._dispatch(arg[1])
            if i < len(t.argnames)-1:
                self._write(', ')
        self._write(")")
        if self._single_func:
            self._do_indent = False
        self._enter()
        self._dispatch(t.code)
        self._leave()
        self._do_indent = True

    def _Getattr(self, t):
        """ Handle getting an attribute of an object
        """
        if isinstance(t.expr, (Div, Mul, Sub, Add)):
            self._write('(')
            self._dispatch(t.expr)
            self._write(')')
        else:
            self._dispatch(t.expr)
            
        self._write('.'+t.attrname)
        
    def _If(self, t):
        self._fill()
        
        for i, (compare,code) in enumerate(t.tests):
            if i == 0:
                self._write("if ")
            else:
                self._write("elif ")
            self._dispatch(compare)
            self._enter()
            self._fill()
            self._dispatch(code)
            self._leave()
            self._write("\n")

        if t.else_ is not None:
            self._write("else")
            self._enter()
            self._fill()
            self._dispatch(t.else_)
            self._leave()
            self._write("\n")
            
    def _IfExp(self, t):
        self._dispatch(t.then)
        self._write(" if ")
        self._dispatch(t.test)

        if t.else_ is not None:
            self._write(" else (")
            self._dispatch(t.else_)
            self._write(")")

    def _Import(self, t):
        """ Handle "import xyz.foo".
        """
        self._fill("import ")
        
        for i, (name,asname) in enumerate(t.names):
            if i != 0:
                self._write(", ")
            self._write(name)
            if asname is not None:
                self._write(" as "+asname)

    def _Keyword(self, t):
        """ Keyword value assignment within function calls and definitions.
        """
        self._write(t.name)
        self._write("=")
        self._dispatch(t.expr)
        
    def _List(self, t):
        self._write("[")
        for  i,node in enumerate(t.nodes):
            self._dispatch(node)
            if i < len(t.nodes)-1:
                self._write(", ")
        self._write("]")

    def _Module(self, t):
        if t.doc is not None:
            self._dispatch(t.doc)
        self._dispatch(t.node)

    def _Mul(self, t):
        self.__binary_op(t, '*')

    def _Name(self, t):
        self._write(t.name)

    def _NoneType(self, t):
        self._write("None")
        
    def _Not(self, t):
        self._write('not (')
        self._dispatch(t.expr)
        self._write(')')
        
    def _Or(self, t):
        self._write(" (")
        for i, node in enumerate(t.nodes):
            self._dispatch(node)
            if i != len(t.nodes)-1:
                self._write(") or (")
        self._write(")")
                
    def _Pass(self, t):
        self._write("pass\n")

    def _Printnl(self, t):
        self._fill("print ")
        if t.dest:
            self._write(">> ")
            self._dispatch(t.dest)
            self._write(", ")
        comma = False
        for node in t.nodes:
            if comma: self._write(', ')
            else: comma = True
            self._dispatch(node)

    def _Power(self, t):
        self.__binary_op(t, '**')

    def _Return(self, t):
        self._fill("return ")
        if t.value:
            if isinstance(t.value, Tuple):
                text = ', '.join([ name.name for name in t.value.asList() ])
                self._write(text)
            else:
                self._dispatch(t.value)
            if not self._do_indent:
                self._write('; ')

    def _Slice(self, t):
        self._dispatch(t.expr)
        self._write("[")
        if t.lower:
            self._dispatch(t.lower)
        self._write(":")
        if t.upper:
            self._dispatch(t.upper)
        #if t.step:
        #    self._write(":")
        #    self._dispatch(t.step)
        self._write("]")

    def _Sliceobj(self, t):
        for i, node in enumerate(t.nodes):
            if i != 0:
                self._write(":")
            if not (isinstance(node, Const) and node.value is None):
                self._dispatch(node)

    def _Stmt(self, tree):
        for node in tree.nodes:
            self._dispatch(node)

    def _Sub(self, t):
        self.__binary_op(t, '-')

    def _Subscript(self, t):
        self._dispatch(t.expr)
        self._write("[")
        for i, value in enumerate(t.subs):
            if i != 0:
                self._write(",")
            self._dispatch(value)
        self._write("]")

    def _TryExcept(self, t):
        self._fill("try")
        self._enter()
        self._dispatch(t.body)
        self._leave()

        for handler in t.handlers:
            self._fill('except ')
            self._dispatch(handler[0])
            if handler[1] is not None:
                self._write(', ')
                self._dispatch(handler[1])
            self._enter()
            self._dispatch(handler[2])
            self._leave()
            
        if t.else_:
            self._fill("else")
            self._enter()
            self._dispatch(t.else_)
            self._leave()

    def _Tuple(self, t):

        if not t.nodes:
            # Empty tuple.
            self._write("()")
        else:
            self._write("(")

            # _write each elements, separated by a comma.
            for element in t.nodes[:-1]:
                self._dispatch(element)
                self._write(", ")

            # Handle the last one without writing comma
            last_element = t.nodes[-1]
            self._dispatch(last_element)

            self._write(")")
            
    def _UnaryAdd(self, t):
        self._write("+")
        self._dispatch(t.expr)
        
    def _UnarySub(self, t):
        self._write("-")
        self._dispatch(t.expr)        

    def _With(self, t):
        self._fill('with ')
        self._dispatch(t.expr)
        if t.vars:
            self._write(' as ')
            self._dispatch(t.vars.name)
        self._enter()
        self._dispatch(t.body)
        self._leave()
        self._write('\n')
        
    def _int(self, t):
        self._write(repr(t))

    def __binary_op(self, t, symbol):
        # Check if parenthesis are needed on left side and then dispatch
        has_paren = False
        left_class = str(t.left.__class__)
        if (left_class in op_precedence.keys() and
            op_precedence[left_class] < op_precedence[str(t.__class__)]):
            has_paren = True
        if has_paren:
            self._write('(')
        self._dispatch(t.left)
        if has_paren:
            self._write(')')
        # Write the appropriate symbol for operator
        self._write(symbol)
        # Check if parenthesis are needed on the right side and then dispatch
        has_paren = False
        right_class = str(t.right.__class__)
        if (right_class in op_precedence.keys() and
            op_precedence[right_class] < op_precedence[str(t.__class__)]):
            has_paren = True
        if has_paren:
            self._write('(')
        self._dispatch(t.right)
        if has_paren:
            self._write(')')

    def _float(self, t):
        # if t is 0.1, str(t)->'0.1' while repr(t)->'0.1000000000001'
        # We prefer str here.
        self._write(str(t))

    def _str(self, t):
        self._write(repr(t))
        
    def _tuple(self, t):
        self._write(str(t))

    #########################################################################
    # These are the methods from the _ast modules unparse.
    #
    # As our needs to handle more advanced code increase, we may want to
    # modify some of the methods below so that they work for compiler.ast.
    #########################################################################

#    # stmt
#    def _Expr(self, tree):
#        self._fill()
#        self._dispatch(tree.value)
#
#    def _Import(self, t):
#        self._fill("import ")
#        first = True
#        for a in t.names:
#            if first:
#                first = False
#            else:
#                self._write(", ")
#            self._write(a.name)
#            if a.asname:
#                self._write(" as "+a.asname)
#
##    def _ImportFrom(self, t):
##        self._fill("from ")
##        self._write(t.module)
##        self._write(" import ")
##        for i, a in enumerate(t.names):
##            if i == 0:
##                self._write(", ")
##            self._write(a.name)
##            if a.asname:
##                self._write(" as "+a.asname)
##        # XXX(jpe) what is level for?
##
#
#    def _Break(self, t):
#        self._fill("break")
#
#    def _Continue(self, t):
#        self._fill("continue")
#
#    def _Delete(self, t):
#        self._fill("del ")
#        self._dispatch(t.targets)
#
#    def _Assert(self, t):
#        self._fill("assert ")
#        self._dispatch(t.test)
#        if t.msg:
#            self._write(", ")
#            self._dispatch(t.msg)
#
#    def _Exec(self, t):
#        self._fill("exec ")
#        self._dispatch(t.body)
#        if t.globals:
#            self._write(" in ")
#            self._dispatch(t.globals)
#        if t.locals:
#            self._write(", ")
#            self._dispatch(t.locals)
#
#    def _Print(self, t):
#        self._fill("print ")
#        do_comma = False
#        if t.dest:
#            self._write(">>")
#            self._dispatch(t.dest)
#            do_comma = True
#        for e in t.values:
#            if do_comma:self._write(", ")
#            else:do_comma=True
#            self._dispatch(e)
#        if not t.nl:
#            self._write(",")
#
#    def _Global(self, t):
#        self._fill("global")
#        for i, n in enumerate(t.names):
#            if i != 0:
#                self._write(",")
#            self._write(" " + n)
#
#    def _Yield(self, t):
#        self._fill("yield")
#        if t.value:
#            self._write(" (")
#            self._dispatch(t.value)
#            self._write(")")
#
#    def _Raise(self, t):
#        self._fill('raise ')
#        if t.type:
#            self._dispatch(t.type)
#        if t.inst:
#            self._write(", ")
#            self._dispatch(t.inst)
#        if t.tback:
#            self._write(", ")
#            self._dispatch(t.tback)
#
#
#    def _TryFinally(self, t):
#        self._fill("try")
#        self._enter()
#        self._dispatch(t.body)
#        self._leave()
#
#        self._fill("finally")
#        self._enter()
#        self._dispatch(t.finalbody)
#        self._leave()
#
#    def _excepthandler(self, t):
#        self._fill("except ")
#        if t.type:
#            self._dispatch(t.type)
#        if t.name:
#            self._write(", ")
#            self._dispatch(t.name)
#        self._enter()
#        self._dispatch(t.body)
#        self._leave()
#
#    def _ClassDef(self, t):
#        self._write("\n")
#        self._fill("class "+t.name)
#        if t.bases:
#            self._write("(")
#            for a in t.bases:
#                self._dispatch(a)
#                self._write(", ")
#            self._write(")")
#        self._enter()
#        self._dispatch(t.body)
#        self._leave()
#
#    def _FunctionDef(self, t):
#        self._write("\n")
#        for deco in t.decorators:
#            self._fill("@")
#            self._dispatch(deco)
#        self._fill("def "+t.name + "(")
#        self._dispatch(t.args)
#        self._write(")")
#        self._enter()
#        self._dispatch(t.body)
#        self._leave()
#
#    def _For(self, t):
#        self._fill("for ")
#        self._dispatch(t.target)
#        self._write(" in ")
#        self._dispatch(t.iter)
#        self._enter()
#        self._dispatch(t.body)
#        self._leave()
#        if t.orelse:
#            self._fill("else")
#            self._enter()
#            self._dispatch(t.orelse)
#            self._leave
#
#    def _While(self, t):
#        self._fill("while ")
#        self._dispatch(t.test)
#        self._enter()
#        self._dispatch(t.body)
#        self._leave()
#        if t.orelse:
#            self._fill("else")
#            self._enter()
#            self._dispatch(t.orelse)
#            self._leave
#
#    # expr
#    def _Str(self, tree):
#        self._write(repr(tree.s))
##
#    def _Repr(self, t):
#        self._write("`")
#        self._dispatch(t.value)
#        self._write("`")
#
#    def _Num(self, t):
#        self._write(repr(t.n))
#
#    def _ListComp(self, t):
#        self._write("[")
#        self._dispatch(t.elt)
#        for gen in t.generators:
#            self._dispatch(gen)
#        self._write("]")
#
#    def _GeneratorExp(self, t):
#        self._write("(")
#        self._dispatch(t.elt)
#        for gen in t.generators:
#            self._dispatch(gen)
#        self._write(")")
#
#    def _comprehension(self, t):
#        self._write(" for ")
#        self._dispatch(t.target)
#        self._write(" in ")
#        self._dispatch(t.iter)
#        for if_clause in t.ifs:
#            self._write(" if ")
#            self._dispatch(if_clause)
#
#    def _IfExp(self, t):
#        self._dispatch(t.body)
#        self._write(" if ")
#        self._dispatch(t.test)
#        if t.orelse:
#            self._write(" else ")
#            self._dispatch(t.orelse)
#
#    unop = {"Invert":"~", "Not": "not", "UAdd":"+", "USub":"-"}
#    def _UnaryOp(self, t):
#        self._write(self.unop[t.op.__class__.__name__])
#        self._write("(")
#        self._dispatch(t.operand)
#        self._write(")")
#
#    binop = { "Add":"+", "Sub":"-", "Mult":"*", "Div":"/", "Mod":"%",
#                    "LShift":">>", "RShift":"<<", "BitOr":"|", "BitXor":"^", "BitAnd":"&",
#                    "FloorDiv":"//", "Pow": "**"}
#    def _BinOp(self, t):
#        self._write("(")
#        self._dispatch(t.left)
#        self._write(")" + self.binop[t.op.__class__.__name__] + "(")
#        self._dispatch(t.right)
#        self._write(")")
#
#    boolops = {_ast.And: 'and', _ast.Or: 'or'}
#    def _BoolOp(self, t):
#        self._write("(")
#        self._dispatch(t.values[0])
#        for v in t.values[1:]:
#            self._write(" %s " % self.boolops[t.op.__class__])
#            self._dispatch(v)
#        self._write(")")
#
#    def _Attribute(self,t):
#        self._dispatch(t.value)
#        self._write(".")
#        self._write(t.attr)
#
##    def _Call(self, t):
##        self._dispatch(t.func)
##        self._write("(")
##        comma = False
##        for e in t.args:
##            if comma: self._write(", ")
##            else: comma = True
##            self._dispatch(e)
##        for e in t.keywords:
##            if comma: self._write(", ")
##            else: comma = True
##            self._dispatch(e)
##        if t.starargs:
##            if comma: self._write(", ")
##            else: comma = True
##            self._write("*")
##            self._dispatch(t.starargs)
##        if t.kwargs:
##            if comma: self._write(", ")
##            else: comma = True
##            self._write("**")
##            self._dispatch(t.kwargs)
##        self._write(")")
#
#    # slice
#    def _Index(self, t):
#        self._dispatch(t.value)
#
#    def _ExtSlice(self, t):
#        for i, d in enumerate(t.dims):
#            if i != 0:
#                self._write(': ')
#            self._dispatch(d)
#
#    # others
#    def _arguments(self, t):
#        first = True
#        nonDef = len(t.args)-len(t.defaults)
#        for a in t.args[0:nonDef]:
#            if first:first = False
#            else: self._write(", ")
#            self._dispatch(a)
#        for a,d in zip(t.args[nonDef:], t.defaults):
#            if first:first = False
#            else: self._write(", ")
#            self._dispatch(a),
#            self._write("=")
#            self._dispatch(d)
#        if t.vararg:
#            if first:first = False
#            else: self._write(", ")
#            self._write("*"+t.vararg)
#        if t.kwarg:
#            if first:first = False
#            else: self._write(", ")
#            self._write("**"+t.kwarg)
#
##    def _keyword(self, t):
##        self._write(t.arg)
##        self._write("=")
##        self._dispatch(t.value)
#
#    def _Lambda(self, t):
#        self._write("lambda ")
#        self._dispatch(t.args)
#        self._write(": ")
#        self._dispatch(t.body)



