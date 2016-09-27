# -*- encoding:utf-8 -*-
from __future__ import division, absolute_import, print_function

import sys, textwrap

from numpydoc.docscrape import NumpyDocString, FunctionDoc, ClassDoc
from numpydoc.docscrape_sphinx import SphinxDocString, SphinxClassDoc
from nose.tools import *

if sys.version_info[0] >= 3:
    sixu = lambda s: s
else:
    sixu = lambda s: unicode(s, 'unicode_escape')


doc_txt = '''\
  numpy.multivariate_normal(mean, cov, shape=None, spam=None)

  Draw values from a multivariate normal distribution with specified
  mean and covariance.

  The multivariate normal or Gaussian distribution is a generalisation
  of the one-dimensional normal distribution to higher dimensions.

  Parameters
  ----------
  mean : (N,) ndarray
      Mean of the N-dimensional distribution.

      .. math::

         (1+2+3)/3

  cov : (N, N) ndarray
      Covariance matrix of the distribution.
  shape : tuple of ints
      Given a shape of, for example, (m,n,k), m*n*k samples are
      generated, and packed in an m-by-n-by-k arrangement.  Because
      each sample is N-dimensional, the output shape is (m,n,k,N).

  Returns
  -------
  out : ndarray
      The drawn samples, arranged according to `shape`.  If the
      shape given is (m,n,...), then the shape of `out` is is
      (m,n,...,N).

      In other words, each entry ``out[i,j,...,:]`` is an N-dimensional
      value drawn from the distribution.
  list of str
      This is not a real return value.  It exists to test
      anonymous return values.

  Other Parameters
  ----------------
  spam : parrot
      A parrot off its mortal coil.

  Raises
  ------
  RuntimeError
      Some error

  Warns
  -----
  RuntimeWarning
      Some warning

  Warnings
  --------
  Certain warnings apply.

  Notes
  -----
  Instead of specifying the full covariance matrix, popular
  approximations include:

    - Spherical covariance (`cov` is a multiple of the identity matrix)
    - Diagonal covariance (`cov` has non-negative elements only on the diagonal)

  This geometrical property can be seen in two dimensions by plotting
  generated data-points:

  >>> mean = [0,0]
  >>> cov = [[1,0],[0,100]] # diagonal covariance, points lie on x or y-axis

  >>> x,y = multivariate_normal(mean,cov,5000).T
  >>> plt.plot(x,y,'x'); plt.axis('equal'); plt.show()

  Note that the covariance matrix must be symmetric and non-negative
  definite.

  References
  ----------
  .. [1] A. Papoulis, "Probability, Random Variables, and Stochastic
         Processes," 3rd ed., McGraw-Hill Companies, 1991
  .. [2] R.O. Duda, P.E. Hart, and D.G. Stork, "Pattern Classification,"
         2nd ed., Wiley, 2001.

  See Also
  --------
  some, other, funcs
  otherfunc : relationship

  Examples
  --------
  >>> mean = (1,2)
  >>> cov = [[1,0],[1,0]]
  >>> x = multivariate_normal(mean,cov,(3,3))
  >>> print x.shape
  (3, 3, 2)

  The following is probably true, given that 0.6 is roughly twice the
  standard deviation:

  >>> print list( (x[0,0,:] - mean) < 0.6 )
  [True, True]

  .. index:: random
     :refguide: random;distributions, random;gauss

  '''
doc = NumpyDocString(doc_txt)

doc_yields_txt = """
Test generator

Yields
------
a : int
    The number of apples.
b : int
    The number of bananas.
int
    The number of unknowns.
"""
doc_yields = NumpyDocString(doc_yields_txt)


def test_signature():
    assert doc['Signature'].startswith('numpy.multivariate_normal(')
    assert doc['Signature'].endswith('spam=None)')

def test_summary():
    assert doc['Summary'][0].startswith('Draw values')
    assert doc['Summary'][-1].endswith('covariance.')

def test_extended_summary():
    assert doc['Extended Summary'][0].startswith('The multivariate normal')

def test_parameters():
    assert_equal(len(doc['Parameters']), 3)
    assert_equal([n for n,_,_ in doc['Parameters']], ['mean','cov','shape'])

    arg, arg_type, desc = doc['Parameters'][1]
    assert_equal(arg_type, '(N, N) ndarray')
    assert desc[0].startswith('Covariance matrix')
    assert doc['Parameters'][0][-1][-2] == '   (1+2+3)/3'

def test_other_parameters():
    assert_equal(len(doc['Other Parameters']), 1)
    assert_equal([n for n,_,_ in doc['Other Parameters']], ['spam'])
    arg, arg_type, desc = doc['Other Parameters'][0]
    assert_equal(arg_type, 'parrot')
    assert desc[0].startswith('A parrot off its mortal coil')

def test_returns():
    assert_equal(len(doc['Returns']), 2)
    arg, arg_type, desc = doc['Returns'][0]
    assert_equal(arg, 'out')
    assert_equal(arg_type, 'ndarray')
    assert desc[0].startswith('The drawn samples')
    assert desc[-1].endswith('distribution.')

    arg, arg_type, desc = doc['Returns'][1]
    assert_equal(arg, 'list of str')
    assert_equal(arg_type, '')
    assert desc[0].startswith('This is not a real')
    assert desc[-1].endswith('anonymous return values.')

def test_yields():
    section = doc_yields['Yields']
    assert_equal(len(section), 3)
    truth = [('a', 'int', 'apples.'),
             ('b', 'int', 'bananas.'),
             ('int', '', 'unknowns.')]
    for (arg, arg_type, desc), (arg_, arg_type_, end) in zip(section, truth):
        assert_equal(arg, arg_)
        assert_equal(arg_type, arg_type_)
        assert desc[0].startswith('The number of')
        assert desc[0].endswith(end)

def test_returnyield():
    doc_text = """
Test having returns and yields.

Returns
-------
int
    The number of apples.

Yields
------
a : int
    The number of apples.
b : int
    The number of bananas.

"""
    assert_raises(ValueError, NumpyDocString, doc_text)

def test_notes():
    assert doc['Notes'][0].startswith('Instead')
    assert doc['Notes'][-1].endswith('definite.')
    assert_equal(len(doc['Notes']), 17)

def test_references():
    assert doc['References'][0].startswith('..')
    assert doc['References'][-1].endswith('2001.')

def test_examples():
    assert doc['Examples'][0].startswith('>>>')
    assert doc['Examples'][-1].endswith('True]')

def test_index():
    assert_equal(doc['index']['default'], 'random')
    assert_equal(len(doc['index']), 2)
    assert_equal(len(doc['index']['refguide']), 2)

def non_blank_line_by_line_compare(a,b):
    a = textwrap.dedent(a)
    b = textwrap.dedent(b)
    a = [l.rstrip() for l in a.split('\n') if l.strip()]
    b = [l.rstrip() for l in b.split('\n') if l.strip()]
    for n,line in enumerate(a):
        if not line == b[n]:
            raise AssertionError("Lines %s of a and b differ: "
                                 "\n>>> %s\n<<< %s\n" %
                                 (n,line,b[n]))
def test_str():
    # doc_txt has the order of Notes and See Also sections flipped.
    # This should be handled automatically, and so, one thing this test does
    # is to make sure that See Also precedes Notes in the output.
    non_blank_line_by_line_compare(str(doc),
"""numpy.multivariate_normal(mean, cov, shape=None, spam=None)

Draw values from a multivariate normal distribution with specified
mean and covariance.

The multivariate normal or Gaussian distribution is a generalisation
of the one-dimensional normal distribution to higher dimensions.

Parameters
----------
mean : (N,) ndarray
    Mean of the N-dimensional distribution.

    .. math::

       (1+2+3)/3

cov : (N, N) ndarray
    Covariance matrix of the distribution.
shape : tuple of ints
    Given a shape of, for example, (m,n,k), m*n*k samples are
    generated, and packed in an m-by-n-by-k arrangement.  Because
    each sample is N-dimensional, the output shape is (m,n,k,N).

Returns
-------
out : ndarray
    The drawn samples, arranged according to `shape`.  If the
    shape given is (m,n,...), then the shape of `out` is is
    (m,n,...,N).

    In other words, each entry ``out[i,j,...,:]`` is an N-dimensional
    value drawn from the distribution.
list of str
    This is not a real return value.  It exists to test
    anonymous return values.

Other Parameters
----------------
spam : parrot
    A parrot off its mortal coil.

Raises
------
RuntimeError
    Some error

Warns
-----
RuntimeWarning
    Some warning

Warnings
--------
Certain warnings apply.

See Also
--------
`some`_, `other`_, `funcs`_

`otherfunc`_
    relationship

Notes
-----
Instead of specifying the full covariance matrix, popular
approximations include:

  - Spherical covariance (`cov` is a multiple of the identity matrix)
  - Diagonal covariance (`cov` has non-negative elements only on the diagonal)

This geometrical property can be seen in two dimensions by plotting
generated data-points:

>>> mean = [0,0]
>>> cov = [[1,0],[0,100]] # diagonal covariance, points lie on x or y-axis

>>> x,y = multivariate_normal(mean,cov,5000).T
>>> plt.plot(x,y,'x'); plt.axis('equal'); plt.show()

Note that the covariance matrix must be symmetric and non-negative
definite.

References
----------
.. [1] A. Papoulis, "Probability, Random Variables, and Stochastic
       Processes," 3rd ed., McGraw-Hill Companies, 1991
.. [2] R.O. Duda, P.E. Hart, and D.G. Stork, "Pattern Classification,"
       2nd ed., Wiley, 2001.

Examples
--------
>>> mean = (1,2)
>>> cov = [[1,0],[1,0]]
>>> x = multivariate_normal(mean,cov,(3,3))
>>> print x.shape
(3, 3, 2)

The following is probably true, given that 0.6 is roughly twice the
standard deviation:

>>> print list( (x[0,0,:] - mean) < 0.6 )
[True, True]

.. index:: random
   :refguide: random;distributions, random;gauss""")


def test_yield_str():
    non_blank_line_by_line_compare(str(doc_yields),
"""Test generator

Yields
------
a : int
    The number of apples.
b : int
    The number of bananas.
int
    The number of unknowns.

.. index:: """)


def test_sphinx_str():
    sphinx_doc = SphinxDocString(doc_txt)
    non_blank_line_by_line_compare(str(sphinx_doc),
"""
.. index:: random
   single: random;distributions, random;gauss

Draw values from a multivariate normal distribution with specified
mean and covariance.

The multivariate normal or Gaussian distribution is a generalisation
of the one-dimensional normal distribution to higher dimensions.

:Parameters:

    **mean** : (N,) ndarray

        Mean of the N-dimensional distribution.

        .. math::

           (1+2+3)/3

    **cov** : (N, N) ndarray

        Covariance matrix of the distribution.

    **shape** : tuple of ints

        Given a shape of, for example, (m,n,k), m*n*k samples are
        generated, and packed in an m-by-n-by-k arrangement.  Because
        each sample is N-dimensional, the output shape is (m,n,k,N).

:Returns:

    **out** : ndarray

        The drawn samples, arranged according to `shape`.  If the
        shape given is (m,n,...), then the shape of `out` is is
        (m,n,...,N).

        In other words, each entry ``out[i,j,...,:]`` is an N-dimensional
        value drawn from the distribution.

    list of str

        This is not a real return value.  It exists to test
        anonymous return values.

:Other Parameters:

    **spam** : parrot

        A parrot off its mortal coil.

:Raises:

    **RuntimeError**

        Some error

:Warns:

    **RuntimeWarning**

        Some warning

.. warning::

    Certain warnings apply.

.. seealso::

    :obj:`some`, :obj:`other`, :obj:`funcs`

    :obj:`otherfunc`
        relationship

.. rubric:: Notes

Instead of specifying the full covariance matrix, popular
approximations include:

  - Spherical covariance (`cov` is a multiple of the identity matrix)
  - Diagonal covariance (`cov` has non-negative elements only on the diagonal)

This geometrical property can be seen in two dimensions by plotting
generated data-points:

>>> mean = [0,0]
>>> cov = [[1,0],[0,100]] # diagonal covariance, points lie on x or y-axis

>>> x,y = multivariate_normal(mean,cov,5000).T
>>> plt.plot(x,y,'x'); plt.axis('equal'); plt.show()

Note that the covariance matrix must be symmetric and non-negative
definite.

.. rubric:: References

.. [1] A. Papoulis, "Probability, Random Variables, and Stochastic
       Processes," 3rd ed., McGraw-Hill Companies, 1991
.. [2] R.O. Duda, P.E. Hart, and D.G. Stork, "Pattern Classification,"
       2nd ed., Wiley, 2001.

.. only:: latex

   [1]_, [2]_

.. rubric:: Examples

>>> mean = (1,2)
>>> cov = [[1,0],[1,0]]
>>> x = multivariate_normal(mean,cov,(3,3))
>>> print x.shape
(3, 3, 2)

The following is probably true, given that 0.6 is roughly twice the
standard deviation:

>>> print list( (x[0,0,:] - mean) < 0.6 )
[True, True]
""")


def test_sphinx_yields_str():
    sphinx_doc = SphinxDocString(doc_yields_txt)
    non_blank_line_by_line_compare(str(sphinx_doc),
"""Test generator

:Yields:

    **a** : int

        The number of apples.

    **b** : int

        The number of bananas.

    int

        The number of unknowns.
""")


doc2 = NumpyDocString("""
    Returns array of indices of the maximum values of along the given axis.

    Parameters
    ----------
    a : {array_like}
        Array to look in.
    axis : {None, integer}
        If None, the index is into the flattened array, otherwise along
        the specified axis""")

def test_parameters_without_extended_description():
    assert_equal(len(doc2['Parameters']), 2)

doc3 = NumpyDocString("""
    my_signature(*params, **kwds)

    Return this and that.
    """)

def test_escape_stars():
    signature = str(doc3).split('\n')[0]
    assert_equal(signature, 'my_signature(\*params, \*\*kwds)')

    def my_func(a, b, **kwargs):
        pass

    fdoc = FunctionDoc(func=my_func)
    assert_equal(fdoc['Signature'], 'my_func(a, b, \*\*kwargs)')

doc4 = NumpyDocString(
    """a.conj()

    Return an array with all complex-valued elements conjugated.""")

def test_empty_extended_summary():
    assert_equal(doc4['Extended Summary'], [])

doc5 = NumpyDocString(
    """
    a.something()

    Raises
    ------
    LinAlgException
        If array is singular.

    Warns
    -----
    SomeWarning
        If needed
    """)

def test_raises():
    assert_equal(len(doc5['Raises']), 1)
    name,_,desc = doc5['Raises'][0]
    assert_equal(name,'LinAlgException')
    assert_equal(desc,['If array is singular.'])

def test_warns():
    assert_equal(len(doc5['Warns']), 1)
    name,_,desc = doc5['Warns'][0]
    assert_equal(name,'SomeWarning')
    assert_equal(desc,['If needed'])

def test_see_also():
    doc6 = NumpyDocString(
    """
    z(x,theta)

    See Also
    --------
    func_a, func_b, func_c
    func_d : some equivalent func
    foo.func_e : some other func over
             multiple lines
    func_f, func_g, :meth:`func_h`, func_j,
    func_k
    :obj:`baz.obj_q`
    :class:`class_j`: fubar
        foobar
    """)

    assert len(doc6['See Also']) == 12
    for func, desc, role in doc6['See Also']:
        if func in ('func_a', 'func_b', 'func_c', 'func_f',
                    'func_g', 'func_h', 'func_j', 'func_k', 'baz.obj_q'):
            assert(not desc)
        else:
            assert(desc)

        if func == 'func_h':
            assert role == 'meth'
        elif func == 'baz.obj_q':
            assert role == 'obj'
        elif func == 'class_j':
            assert role == 'class'
        else:
            assert role is None

        if func == 'func_d':
            assert desc == ['some equivalent func']
        elif func == 'foo.func_e':
            assert desc == ['some other func over', 'multiple lines']
        elif func == 'class_j':
            assert desc == ['fubar', 'foobar']

def test_see_also_print():
    class Dummy(object):
        """
        See Also
        --------
        func_a, func_b
        func_c : some relationship
                 goes here
        func_d
        """
        pass

    obj = Dummy()
    s = str(FunctionDoc(obj, role='func'))
    assert(':func:`func_a`, :func:`func_b`' in s)
    assert('    some relationship' in s)
    assert(':func:`func_d`' in s)

doc7 = NumpyDocString("""

        Doc starts on second line.

        """)

def test_empty_first_line():
    assert doc7['Summary'][0].startswith('Doc starts')


def test_no_summary():
    str(SphinxDocString("""
    Parameters
    ----------"""))


def test_unicode():
    doc = SphinxDocString("""
    öäöäöäöäöåååå

    öäöäöäööäååå

    Parameters
    ----------
    ååå : äää
        ööö

    Returns
    -------
    ååå : ööö
        äää

    """)
    assert isinstance(doc['Summary'][0], str)
    assert doc['Summary'][0] == 'öäöäöäöäöåååå'

def test_plot_examples():
    cfg = dict(use_plots=True)

    doc = SphinxDocString("""
    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> plt.plot([1,2,3],[4,5,6])
    >>> plt.show()
    """, config=cfg)
    assert 'plot::' in str(doc), str(doc)

    doc = SphinxDocString("""
    Examples
    --------
    .. plot::

       import matplotlib.pyplot as plt
       plt.plot([1,2,3],[4,5,6])
       plt.show()
    """, config=cfg)
    assert str(doc).count('plot::') == 1, str(doc)

def test_class_members():

    class Dummy(object):
        """
        Dummy class.

        """
        def spam(self, a, b):
            """Spam\n\nSpam spam."""
            pass
        def ham(self, c, d):
            """Cheese\n\nNo cheese."""
            pass
        @property
        def spammity(self):
            """Spammity index"""
            return 0.95

        class Ignorable(object):
            """local class, to be ignored"""
            pass

    for cls in (ClassDoc, SphinxClassDoc):
        doc = cls(Dummy, config=dict(show_class_members=False))
        assert 'Methods' not in str(doc), (cls, str(doc))
        assert 'spam' not in str(doc), (cls, str(doc))
        assert 'ham' not in str(doc), (cls, str(doc))
        assert 'spammity' not in str(doc), (cls, str(doc))
        assert 'Spammity index' not in str(doc), (cls, str(doc))

        doc = cls(Dummy, config=dict(show_class_members=True))
        assert 'Methods' in str(doc), (cls, str(doc))
        assert 'spam' in str(doc), (cls, str(doc))
        assert 'ham' in str(doc), (cls, str(doc))
        assert 'spammity' in str(doc), (cls, str(doc))

        if cls is SphinxClassDoc:
            assert '.. autosummary::' in str(doc), str(doc)
        else:
            assert 'Spammity index' in str(doc), str(doc)

    class SubDummy(Dummy):
        """
        Subclass of Dummy class.

        """
        def ham(self, c, d):
            """Cheese\n\nNo cheese.\nOverloaded Dummy.ham"""
            pass

        def bar(self, a, b):
            """Bar\n\nNo bar"""
            pass

    for cls in (ClassDoc, SphinxClassDoc):
        doc = cls(SubDummy, config=dict(show_class_members=True,
                                        show_inherited_class_members=False))
        assert 'Methods' in str(doc), (cls, str(doc))
        assert 'spam' not in str(doc), (cls, str(doc))
        assert 'ham' in str(doc), (cls, str(doc))
        assert 'bar' in str(doc), (cls, str(doc))
        assert 'spammity' not in str(doc), (cls, str(doc))

        if cls is SphinxClassDoc:
            assert '.. autosummary::' in str(doc), str(doc)
        else:
            assert 'Spammity index' not in str(doc), str(doc)

        doc = cls(SubDummy, config=dict(show_class_members=True,
                                        show_inherited_class_members=True))
        assert 'Methods' in str(doc), (cls, str(doc))
        assert 'spam' in str(doc), (cls, str(doc))
        assert 'ham' in str(doc), (cls, str(doc))
        assert 'bar' in str(doc), (cls, str(doc))
        assert 'spammity' in str(doc), (cls, str(doc))

        if cls is SphinxClassDoc:
            assert '.. autosummary::' in str(doc), str(doc)
        else:
            assert 'Spammity index' in str(doc), str(doc)

def test_duplicate_signature():
    # Duplicate function signatures occur e.g. in ufuncs, when the
    # automatic mechanism adds one, and a more detailed comes from the
    # docstring itself.

    doc = NumpyDocString(
    """
    z(x1, x2)

    z(a, theta)
    """)

    assert doc['Signature'].strip() == 'z(a, theta)'


class_doc_txt = """
    Foo

    Parameters
    ----------
    f : callable ``f(t, y, *f_args)``
        Aaa.
    jac : callable ``jac(t, y, *jac_args)``
        Bbb.

    Attributes
    ----------
    t : float
        Current time.
    y : ndarray
        Current variable values.
    x : float
        Some parameter

    Methods
    -------
    a
    b
    c

    Examples
    --------
    For usage examples, see `ode`.
"""

def test_class_members_doc():
    doc = ClassDoc(None, class_doc_txt)
    non_blank_line_by_line_compare(str(doc),
    """
    Foo

    Parameters
    ----------
    f : callable ``f(t, y, *f_args)``
        Aaa.
    jac : callable ``jac(t, y, *jac_args)``
        Bbb.

    Examples
    --------
    For usage examples, see `ode`.

    Attributes
    ----------
    t : float
        Current time.
    y : ndarray
        Current variable values.
    x : float
        Some parameter

    Methods
    -------
    a

    b

    c

    .. index::

    """)

def test_class_members_doc_sphinx():
    class Foo:
        @property
        def x(self):
            """Test attribute"""
            return None

    doc = SphinxClassDoc(Foo, class_doc_txt)
    non_blank_line_by_line_compare(str(doc),
    """
    Foo

    :Parameters:

        **f** : callable ``f(t, y, *f_args)``

            Aaa.

        **jac** : callable ``jac(t, y, *jac_args)``

            Bbb.

    .. rubric:: Examples

    For usage examples, see `ode`.

    .. rubric:: Attributes

    .. autosummary::
       :toctree:

       x

    ===  ==========
      t  (float) Current time.
      y  (ndarray) Current variable values.
    ===  ==========

    .. rubric:: Methods

    ===  ==========
      a
      b
      c
    ===  ==========

    """)

if __name__ == "__main__":
    import nose
    nose.run()
