import pytest
import numpy as np
from numpy.testing import assert_array_equal, assert_equal
from numpy.f2py.crackfortran import markinnerspaces
from . import util
from numpy.f2py import crackfortran
import textwrap


class TestNoSpace(util.F2PyTest):
    # issue gh-15035: add handling for endsubroutine, endfunction with no space
    # between "end" and the block name
    code = """
        subroutine subb(k)
          real(8), intent(inout) :: k(:)
          k=k+1
        endsubroutine

        subroutine subc(w,k)
          real(8), intent(in) :: w(:)
          real(8), intent(out) :: k(size(w))
          k=w+1
        endsubroutine

        function t0(value)
          character value
          character t0
          t0 = value
        endfunction
    """

    def test_module(self):
        k = np.array([1, 2, 3], dtype=np.float64)
        w = np.array([1, 2, 3], dtype=np.float64)
        self.module.subb(k)
        assert_array_equal(k, w + 1)
        self.module.subc([w, k])
        assert_array_equal(k, w + 1)
        assert self.module.t0(23) == b'2'


class TestPublicPrivate():

    def test_defaultPrivate(self, tmp_path):
        f_path = tmp_path / "mod.f90"
        with f_path.open('w') as ff:
            ff.write(textwrap.dedent("""\
            module foo
              private
              integer :: a
              public :: setA
              integer :: b
            contains
              subroutine setA(v)
                integer, intent(in) :: v
                a = v
              end subroutine setA
            end module foo
            """))
        mod = crackfortran.crackfortran([str(f_path)])
        assert len(mod) == 1
        mod = mod[0]
        assert 'private' in mod['vars']['a']['attrspec']
        assert 'public' not in mod['vars']['a']['attrspec']
        assert 'private' in mod['vars']['b']['attrspec']
        assert 'public' not in mod['vars']['b']['attrspec']
        assert 'private' not in mod['vars']['seta']['attrspec']
        assert 'public' in mod['vars']['seta']['attrspec']

    def test_defaultPublic(self, tmp_path):
        f_path = tmp_path / "mod.f90"
        with f_path.open('w') as ff:
            ff.write(textwrap.dedent("""\
            module foo
              public
              integer, private :: a
              public :: setA
            contains
              subroutine setA(v)
                integer, intent(in) :: v
                a = v
              end subroutine setA
            end module foo
            """))
        mod = crackfortran.crackfortran([str(f_path)])
        assert len(mod) == 1
        mod = mod[0]
        assert 'private' in mod['vars']['a']['attrspec']
        assert 'public' not in mod['vars']['a']['attrspec']
        assert 'private' not in mod['vars']['seta']['attrspec']
        assert 'public' in mod['vars']['seta']['attrspec']


class TestExternal(util.F2PyTest):
    # issue gh-17859: add external attribute support
    code = """
        integer(8) function external_as_statement(fcn)
        implicit none
        external fcn
        integer(8) :: fcn
        external_as_statement = fcn(0)
        end

        integer(8) function external_as_attribute(fcn)
        implicit none
        integer(8), external :: fcn
        external_as_attribute = fcn(0)
        end
    """

    def test_external_as_statement(self):
        def incr(x):
            return x + 123
        r = self.module.external_as_statement(incr)
        assert r == 123

    def test_external_as_attribute(self):
        def incr(x):
            return x + 123
        r = self.module.external_as_attribute(incr)
        assert r == 123


class TestCrackFortran(util.F2PyTest):

    suffix = '.f90'

    code = textwrap.dedent("""
      subroutine gh2848( &
        ! first 2 parameters
        par1, par2,&
        ! last 2 parameters
        par3, par4)

        integer, intent(in)  :: par1, par2
        integer, intent(out) :: par3, par4

        par3 = par1
        par4 = par2

      end subroutine gh2848
    """)

    def test_gh2848(self):
        r = self.module.gh2848(1, 2)
        assert r == (1, 2)


class TestMarkinnerspaces():
    # issue #14118: markinnerspaces does not handle multiple quotations

    def test_do_not_touch_normal_spaces(self):
        test_list = ["a ", " a", "a b c", "'abcdefghij'"]
        for i in test_list:
            assert_equal(markinnerspaces(i), i)

    def test_one_relevant_space(self):
        assert_equal(markinnerspaces("a 'b c' \\\' \\\'"), "a 'b@_@c' \\' \\'")
        assert_equal(markinnerspaces(r'a "b c" \" \"'), r'a "b@_@c" \" \"')

    def test_ignore_inner_quotes(self):
        assert_equal(markinnerspaces('a \'b c" " d\' e'),
                     "a 'b@_@c\"@_@\"@_@d' e")
        assert_equal(markinnerspaces('a "b c\' \' d" e'),
                     "a \"b@_@c'@_@'@_@d\" e")

    def test_multiple_relevant_spaces(self):
        assert_equal(markinnerspaces("a 'b c' 'd e'"), "a 'b@_@c' 'd@_@e'")
        assert_equal(markinnerspaces(r'a "b c" "d e"'), r'a "b@_@c" "d@_@e"')


class TestDimSpec(util.F2PyTest):
    """This test suite tests various expressions that are used as dimension
    specifications.

    There exists two usage cases where analyzing dimensions
    specifications are important.

    In the first case, the size of output arrays must be defined based
    on the inputs to a Fortran function. Because Fortran supports
    arbitrary bases for indexing, for instance, `arr(lower:upper)`,
    f2py has to evaluate an expression `upper - lower + 1` where
    `lower` and `upper` are arbitrary expressions of input parameters.
    The evaluation is performed in C, so f2py has to translate Fortran
    expressions to valid C expressions (an alternative approach is
    that a developer specifies the corresponding C expressions in a
    .pyf file).

    In the second case, when user provides an input array with a given
    size but some hidden parameters used in dimensions specifications
    need to be determined based on the input array size. This is a
    harder problem because f2py has to solve the inverse problem: find
    a parameter `p` such that `upper(p) - lower(p) + 1` equals to the
    size of input array. In the case when this equation cannot be
    solved (e.g. because the input array size is wrong), raise an
    error before calling the Fortran function (that otherwise would
    likely crash Python process when the size of input arrays is
    wrong). f2py currently supports this case only when the equation
    is linear with respect to unknown parameter.

    """

    suffix = '.f90'

    code_template = textwrap.dedent("""
      function get_arr_size_{count}(a, n) result (length)
        integer, intent(in) :: n
        integer, dimension({dimspec}), intent(out) :: a
        integer length
        length = size(a)
      end function

      subroutine get_inv_arr_size_{count}(a, n)
        integer :: n
        ! the value of n is computed in f2py wrapper
        !f2py intent(out) n
        integer, dimension({dimspec}), intent(in) :: a
        if (a({first}).gt.0) then
          print*, "a=", a
        endif
      end subroutine
    """)

    linear_dimspecs = [
        "n", "2*n", "2:n", "n/2", "5 - n/2", "3*n:20", "n*(n+1):n*(n+5)",
        "2*n, n"
    ]
    nonlinear_dimspecs = ["2*n:3*n*n+2*n"]
    all_dimspecs = linear_dimspecs + nonlinear_dimspecs

    code = ''
    for count, dimspec in enumerate(all_dimspecs):
        lst = [(d.split(":")[0] if ":" in d else "1") for d in dimspec.split(',')]
        code += code_template.format(
            count=count,
            dimspec=dimspec,
            first=", ".join(lst),
        )

    @pytest.mark.parametrize('dimspec', all_dimspecs)
    def test_array_size(self, dimspec):

        count = self.all_dimspecs.index(dimspec)
        get_arr_size = getattr(self.module, f'get_arr_size_{count}')

        for n in [1, 2, 3, 4, 5]:
            sz, a = get_arr_size(n)
            assert a.size == sz

    @pytest.mark.parametrize('dimspec', all_dimspecs)
    def test_inv_array_size(self, dimspec):

        count = self.all_dimspecs.index(dimspec)
        get_arr_size = getattr(self.module, f'get_arr_size_{count}')
        get_inv_arr_size = getattr(self.module, f'get_inv_arr_size_{count}')

        for n in [1, 2, 3, 4, 5]:
            sz, a = get_arr_size(n)
            if dimspec in self.nonlinear_dimspecs:
                # one must specify n as input, the call we'll ensure
                # that a and n are compatible:
                n1 = get_inv_arr_size(a, n)
            else:
                # in case of linear dependence, n can be determined
                # from the shape of a:
                n1 = get_inv_arr_size(a)
            # n1 may be different from n (for instance, when `a` size
            # is a function of some `n` fraction) but it must produce
            # the same sized array
            sz1, _ = get_arr_size(n1)
            assert sz == sz1, (n, n1, sz, sz1)


class TestModuleDeclaration():
    def test_dependencies(self, tmp_path):
        f_path = tmp_path / "mod.f90"
        with f_path.open('w') as ff:
            ff.write(textwrap.dedent("""\
            module foo
              type bar
                character(len = 4) :: text
              end type bar
              type(bar), parameter :: abar = bar('abar')
            end module foo
            """))
        mod = crackfortran.crackfortran([str(f_path)])
        assert len(mod) == 1
        assert mod[0]['vars']['abar']['='] == "bar('abar')"
