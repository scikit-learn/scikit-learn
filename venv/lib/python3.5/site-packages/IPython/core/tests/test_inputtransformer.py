import tokenize
import nose.tools as nt

from IPython.testing import tools as tt
from IPython.utils import py3compat
u_fmt = py3compat.u_format

from IPython.core import inputtransformer as ipt

def transform_and_reset(transformer):
    transformer = transformer()
    def transform(inp):
        try:
            return transformer.push(inp)
        finally:
            transformer.reset()
    
    return transform

# Transformer tests
def transform_checker(tests, transformer, **kwargs):
    """Utility to loop over test inputs"""
    transformer = transformer(**kwargs)
    try:
        for inp, tr in tests:
            if inp is None:
                out = transformer.reset()
            else:
                out = transformer.push(inp)
            nt.assert_equal(out, tr)
    finally:
        transformer.reset()

# Data for all the syntax tests in the form of lists of pairs of
# raw/transformed input.  We store it here as a global dict so that we can use
# it both within single-function tests and also to validate the behavior of the
# larger objects

syntax = \
  dict(assign_system =
       [(i,py3compat.u_format(o)) for i,o in \
       [(u'a =! ls', "a = get_ipython().getoutput('ls')"),
        (u'b = !ls', "b = get_ipython().getoutput('ls')"),
        (u'c= !ls', "c = get_ipython().getoutput('ls')"),
        (u'd == !ls', u'd == !ls'), # Invalid syntax, but we leave == alone.
        ('x=1', 'x=1'), # normal input is unmodified
        ('    ','    '),  # blank lines are kept intact
        # Tuple unpacking
        (u"a, b = !echo 'a\\nb'", u"a, b = get_ipython().getoutput(\"echo 'a\\\\nb'\")"),
        (u"a,= !echo 'a'", u"a, = get_ipython().getoutput(\"echo 'a'\")"),
        (u"a, *bc = !echo 'a\\nb\\nc'", u"a, *bc = get_ipython().getoutput(\"echo 'a\\\\nb\\\\nc'\")"),
        # Tuple unpacking with regular Python expressions, not our syntax.
        (u"a, b = range(2)", u"a, b = range(2)"),
        (u"a, = range(1)", u"a, = range(1)"),
        (u"a, *bc = range(3)", u"a, *bc = range(3)"),
        ]],

       assign_magic =
       [(i,py3compat.u_format(o)) for i,o in \
       [(u'a =% who', "a = get_ipython().run_line_magic('who', '')"),
        (u'b = %who', "b = get_ipython().run_line_magic('who', '')"),
        (u'c= %ls', "c = get_ipython().run_line_magic('ls', '')"),
        (u'd == %ls', u'd == %ls'), # Invalid syntax, but we leave == alone.
        ('x=1', 'x=1'), # normal input is unmodified
        ('    ','    '),  # blank lines are kept intact
        (u"a, b = %foo", u"a, b = get_ipython().run_line_magic('foo', '')"),
        ]],

       classic_prompt =
       [('>>> x=1', 'x=1'),
        ('x=1', 'x=1'), # normal input is unmodified
        ('    ', '    '),  # blank lines are kept intact
        ],

       ipy_prompt =
       [('In [1]: x=1', 'x=1'),
        ('x=1', 'x=1'), # normal input is unmodified
        ('    ','    '),  # blank lines are kept intact
        ],

       # Tests for the escape transformer to leave normal code alone
       escaped_noesc =
       [ ('    ', '    '),
         ('x=1', 'x=1'),
         ],

       # System calls
       escaped_shell =
       [(i,py3compat.u_format(o)) for i,o in \
       [ (u'!ls', "get_ipython().system('ls')"),
         # Double-escape shell, this means to capture the output of the
         # subprocess and return it
         (u'!!ls', "get_ipython().getoutput('ls')"),
         ]],

       # Help/object info
       escaped_help =
       [(i,py3compat.u_format(o)) for i,o in \
       [ (u'?', 'get_ipython().show_usage()'),
         (u'?x1', "get_ipython().run_line_magic('pinfo', 'x1')"),
         (u'??x2', "get_ipython().run_line_magic('pinfo2', 'x2')"),
         (u'?a.*s', "get_ipython().run_line_magic('psearch', 'a.*s')"),
         (u'?%hist1', "get_ipython().run_line_magic('pinfo', '%hist1')"),
         (u'?%%hist2', "get_ipython().run_line_magic('pinfo', '%%hist2')"),
         (u'?abc = qwe', "get_ipython().run_line_magic('pinfo', 'abc')"),
         ]],

      end_help =
      [(i,py3compat.u_format(o)) for i,o in \
      [ (u'x3?', "get_ipython().run_line_magic('pinfo', 'x3')"),
        (u'x4??', "get_ipython().run_line_magic('pinfo2', 'x4')"),
        (u'%hist1?', "get_ipython().run_line_magic('pinfo', '%hist1')"),
        (u'%hist2??', "get_ipython().run_line_magic('pinfo2', '%hist2')"),
        (u'%%hist3?', "get_ipython().run_line_magic('pinfo', '%%hist3')"),
        (u'%%hist4??', "get_ipython().run_line_magic('pinfo2', '%%hist4')"),
        (u'f*?', "get_ipython().run_line_magic('psearch', 'f*')"),
        (u'ax.*aspe*?', "get_ipython().run_line_magic('psearch', 'ax.*aspe*')"),
        (u'a = abc?', "get_ipython().set_next_input('a = abc');"
                      "get_ipython().run_line_magic('pinfo', 'abc')"),
        (u'a = abc.qe??', "get_ipython().set_next_input('a = abc.qe');"
                          "get_ipython().run_line_magic('pinfo2', 'abc.qe')"),
        (u'a = *.items?', "get_ipython().set_next_input('a = *.items');"
                          "get_ipython().run_line_magic('psearch', '*.items')"),
        (u'plot(a?', "get_ipython().set_next_input('plot(a');"
                     "get_ipython().run_line_magic('pinfo', 'a')"),
        (u'a*2 #comment?', 'a*2 #comment?'),
        ]],

       # Explicit magic calls
       escaped_magic =
       [(i,py3compat.u_format(o)) for i,o in \
       [ (u'%cd', "get_ipython().run_line_magic('cd', '')"),
         (u'%cd /home', "get_ipython().run_line_magic('cd', '/home')"),
         # Backslashes need to be escaped.
         (u'%cd C:\\User', "get_ipython().run_line_magic('cd', 'C:\\\\User')"),
         (u'    %magic', "    get_ipython().run_line_magic('magic', '')"),
         ]],

       # Quoting with separate arguments
       escaped_quote =
       [ (',f', 'f("")'),
         (',f x', 'f("x")'),
         ('  ,f y', '  f("y")'),
         (',f a b', 'f("a", "b")'),
         ],

       # Quoting with single argument
       escaped_quote2 =
       [ (';f', 'f("")'),
         (';f x', 'f("x")'),
         ('  ;f y', '  f("y")'),
         (';f a b', 'f("a b")'),
         ],

       # Simply apply parens
       escaped_paren =
       [ ('/f', 'f()'),
         ('/f x', 'f(x)'),
         ('  /f y', '  f(y)'),
         ('/f a b', 'f(a, b)'),
         ],

       # Check that we transform prompts before other transforms
       mixed =
       [(i,py3compat.u_format(o)) for i,o in \
       [ (u'In [1]: %lsmagic', "get_ipython().run_line_magic('lsmagic', '')"),
         (u'>>> %lsmagic', "get_ipython().run_line_magic('lsmagic', '')"),
         (u'In [2]: !ls', "get_ipython().system('ls')"),
         (u'In [3]: abs?', "get_ipython().run_line_magic('pinfo', 'abs')"),
         (u'In [4]: b = %who', "b = get_ipython().run_line_magic('who', '')"),
         ]],
       )

# multiline syntax examples.  Each of these should be a list of lists, with
# each entry itself having pairs of raw/transformed input.  The union (with
# '\n'.join() of the transformed inputs is what the splitter should produce
# when fed the raw lines one at a time via push.
syntax_ml = \
  dict(classic_prompt =
       [ [('>>> for i in range(10):','for i in range(10):'),
          ('...     print i','    print i'),
          ('... ', ''),
          ],
         [('>>> a="""','a="""'),
          ('... 123"""','123"""'),
          ],
         [('a="""','a="""'),
          ('... 123','123'),
          ('... 456"""','456"""'),
          ],
         [('a="""','a="""'),
          ('>>> 123','123'),
          ('... 456"""','456"""'),
          ],
         [('a="""','a="""'),
          ('123','123'),
          ('... 456"""','... 456"""'),
          ],
         [('....__class__','....__class__'),
         ],
         [('a=5', 'a=5'),
          ('...', ''),
         ],
         [('>>> def f(x):', 'def f(x):'),
          ('...', ''),
          ('...     return x', '    return x'),
          ],
         [('board = """....', 'board = """....'),
          ('....', '....'),
          ('...."""', '...."""'),
          ],
        ],

       ipy_prompt =
       [ [('In [24]: for i in range(10):','for i in range(10):'),
          ('   ....:     print i','    print i'),
          ('   ....: ', ''),
          ],
         [('In [24]: for i in range(10):','for i in range(10):'),
          # Qt console prompts expand with spaces, not dots
          ('    ...:     print i','    print i'),
          ('    ...: ', ''),
          ],
         [('In [24]: for i in range(10):','for i in range(10):'),
          # Sometimes whitespace preceding '...' has been removed
          ('...:     print i','    print i'),
          ('...: ', ''),
          ],
         [('In [24]: for i in range(10):','for i in range(10):'),
          # Space after last continuation prompt has been removed (issue #6674)
          ('...:     print i','    print i'),
          ('...:', ''),
          ],
         [('In [2]: a="""','a="""'),
          ('   ...: 123"""','123"""'),
          ],
         [('a="""','a="""'),
          ('   ...: 123','123'),
          ('   ...: 456"""','456"""'),
          ],
         [('a="""','a="""'),
          ('In [1]: 123','123'),
          ('   ...: 456"""','456"""'),
          ],
         [('a="""','a="""'),
          ('123','123'),
          ('   ...: 456"""','   ...: 456"""'),
          ],
         ],

       multiline_datastructure_prompt =
       [ [('>>> a = [1,','a = [1,'),
          ('... 2]','2]'),
         ],
       ],
        
       multiline_datastructure =
       [ [('b = ("%s"', None),
          ('# comment', None),
          ('%foo )', 'b = ("%s"\n# comment\n%foo )'),
         ],
       ],
       
       multiline_string =
       [ [("'''foo?", None),
          ("bar'''", "'''foo?\nbar'''"),
         ],
       ],
       
       leading_indent =
       [ [('    print "hi"','print "hi"'),
          ],
         [('  for a in range(5):','for a in range(5):'),
          ('    a*2','  a*2'),
          ],
         [('    a="""','a="""'),
          ('    123"""','123"""'),
           ],
         [('a="""','a="""'),
          ('    123"""','    123"""'),
          ],
       ],
       
       cellmagic =
       [ [(u'%%foo a', None),
          (None, u_fmt("get_ipython().run_cell_magic('foo', 'a', '')")),
          ],
         [(u'%%bar 123', None),
          (u'hello', None),
          (None , u_fmt("get_ipython().run_cell_magic('bar', '123', 'hello')")),
          ],
         [(u'a=5', 'a=5'),
          (u'%%cellmagic', '%%cellmagic'),
          ],
       ],
       
       escaped =
       [ [('%abc def \\', None),
          ('ghi', u_fmt("get_ipython().run_line_magic('abc', 'def ghi')")),
          ],
         [('%abc def \\', None),
          ('ghi\\', None),
          (None, u_fmt("get_ipython().run_line_magic('abc', 'def ghi')")),
          ],
       ],
       
       assign_magic =
       [ [(u'a = %bc de \\', None),
          (u'fg', u_fmt("a = get_ipython().run_line_magic('bc', 'de fg')")),
          ],
         [(u'a = %bc de \\', None),
          (u'fg\\', None),
          (None, u_fmt("a = get_ipython().run_line_magic('bc', 'de fg')")),
          ],
       ],
       
       assign_system =
       [ [(u'a = !bc de \\', None),
          (u'fg', u_fmt("a = get_ipython().getoutput('bc de fg')")),
          ],
         [(u'a = !bc de \\', None),
          (u'fg\\', None),
          (None, u_fmt("a = get_ipython().getoutput('bc de fg')")),
          ],
       ],
       )


def test_assign_system():
    tt.check_pairs(transform_and_reset(ipt.assign_from_system), syntax['assign_system'])

def test_assign_magic():
    tt.check_pairs(transform_and_reset(ipt.assign_from_magic), syntax['assign_magic'])

def test_classic_prompt():
    tt.check_pairs(transform_and_reset(ipt.classic_prompt), syntax['classic_prompt'])
    for example in syntax_ml['classic_prompt']:
        transform_checker(example, ipt.classic_prompt)
    for example in syntax_ml['multiline_datastructure_prompt']:
        transform_checker(example, ipt.classic_prompt)

    # Check that we don't transform the second line if the first is obviously
    # IPython syntax
    transform_checker([
        (u'%foo', '%foo'),
        (u'>>> bar', '>>> bar'),
    ], ipt.classic_prompt)


def test_ipy_prompt():
    tt.check_pairs(transform_and_reset(ipt.ipy_prompt), syntax['ipy_prompt'])
    for example in syntax_ml['ipy_prompt']:
        transform_checker(example, ipt.ipy_prompt)

    # Check that we don't transform the second line if we're inside a cell magic
    transform_checker([
        (u'%%foo', '%%foo'),
        (u'In [1]: bar', 'In [1]: bar'),
    ], ipt.ipy_prompt)

def test_assemble_logical_lines():
    tests = \
    [ [(u"a = \\", None),
       (u"123", u"a = 123"),
      ],
      [(u"a = \\", None),  # Test resetting when within a multi-line string
       (u"12 *\\", None),
       (None, u"a = 12 *"),
      ],
      [(u"# foo\\", u"# foo\\"), # Comments can't be continued like this
      ],
    ]
    for example in tests:
        transform_checker(example, ipt.assemble_logical_lines)

def test_assemble_python_lines():
    tests = \
    [ [(u"a = '''", None),
       (u"abc'''", u"a = '''\nabc'''"),
      ],
      [(u"a = '''", None),  # Test resetting when within a multi-line string
       (u"def", None),
       (None, u"a = '''\ndef"),
      ],
      [(u"a = [1,", None),
       (u"2]", u"a = [1,\n2]"),
      ],
      [(u"a = [1,", None),  # Test resetting when within a multi-line string
       (u"2,", None),
       (None, u"a = [1,\n2,"),
      ],
      [(u"a = '''", None),  # Test line continuation within a multi-line string
       (u"abc\\", None),
       (u"def", None),
       (u"'''", u"a = '''\nabc\\\ndef\n'''"),
      ],
    ] + syntax_ml['multiline_datastructure']
    for example in tests:
        transform_checker(example, ipt.assemble_python_lines)


def test_help_end():
    tt.check_pairs(transform_and_reset(ipt.help_end), syntax['end_help'])

def test_escaped_noesc():
    tt.check_pairs(transform_and_reset(ipt.escaped_commands), syntax['escaped_noesc'])


def test_escaped_shell():
    tt.check_pairs(transform_and_reset(ipt.escaped_commands), syntax['escaped_shell'])


def test_escaped_help():
    tt.check_pairs(transform_and_reset(ipt.escaped_commands), syntax['escaped_help'])


def test_escaped_magic():
    tt.check_pairs(transform_and_reset(ipt.escaped_commands), syntax['escaped_magic'])


def test_escaped_quote():
    tt.check_pairs(transform_and_reset(ipt.escaped_commands), syntax['escaped_quote'])


def test_escaped_quote2():
    tt.check_pairs(transform_and_reset(ipt.escaped_commands), syntax['escaped_quote2'])


def test_escaped_paren():
    tt.check_pairs(transform_and_reset(ipt.escaped_commands), syntax['escaped_paren'])


def test_cellmagic():
    for example in syntax_ml['cellmagic']:
        transform_checker(example, ipt.cellmagic)
    
    line_example = [(u'%%bar 123', None),
                    (u'hello', None),
                    (u'' , u_fmt("get_ipython().run_cell_magic('bar', '123', 'hello')")),
                   ]
    transform_checker(line_example, ipt.cellmagic, end_on_blank_line=True)

def test_has_comment():
    tests = [('text', False),
             ('text #comment', True),
             ('text #comment\n', True),
             ('#comment', True),
             ('#comment\n', True),
             ('a = "#string"', False),
             ('a = "#string" # comment', True),
             ('a #comment not "string"', True),
             ]
    tt.check_pairs(ipt.has_comment, tests)

@ipt.TokenInputTransformer.wrap
def decistmt(tokens):
    """Substitute Decimals for floats in a string of statements.

    Based on an example from the tokenize module docs.
    """
    result = []
    for toknum, tokval, _, _, _  in tokens:
        if toknum == tokenize.NUMBER and '.' in tokval:  # replace NUMBER tokens
            for newtok in [
                (tokenize.NAME, 'Decimal'),
                (tokenize.OP, '('),
                (tokenize.STRING, repr(tokval)),
                (tokenize.OP, ')')
            ]:
                yield newtok
        else:
            yield (toknum, tokval)



def test_token_input_transformer():
    tests = [(u'1.2', u_fmt(u"Decimal ('1.2')")),
             (u'"1.2"', u'"1.2"'),
             ]
    tt.check_pairs(transform_and_reset(decistmt), tests)
    ml_tests = \
    [ [(u"a = 1.2; b = '''x", None),
       (u"y'''", u_fmt(u"a =Decimal ('1.2');b ='''x\ny'''")),
      ],
      [(u"a = [1.2,", None),
       (u"3]", u_fmt(u"a =[Decimal ('1.2'),\n3 ]")),
      ],
      [(u"a = '''foo", None),  # Test resetting when within a multi-line string
       (u"bar", None),
       (None, u"a = '''foo\nbar"),
      ],
    ]
    for example in ml_tests:
        transform_checker(example, decistmt)
