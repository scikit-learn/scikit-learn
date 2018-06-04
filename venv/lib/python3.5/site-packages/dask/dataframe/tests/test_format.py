# coding: utf-8
import pandas as pd
from textwrap import dedent

import dask.dataframe as dd
from dask.dataframe.utils import PANDAS_VERSION

if PANDAS_VERSION >= '0.21.0':
    style = """<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
"""
elif PANDAS_VERSION >= '0.20.0':
    style = """<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
"""
else:
    style = ""


def test_repr():
    df = pd.DataFrame({'x': list(range(100))})
    ddf = dd.from_pandas(df, 3)

    for x in [ddf, ddf.index, ddf.x]:
        assert type(x).__name__ in repr(x)
        assert str(x.npartitions) in repr(x)


def test_repr_meta_mutation():
    # Check that the repr changes when meta changes
    df = pd.DataFrame({'a': range(5),
                       'b': ['a', 'b', 'c', 'd', 'e']})
    ddf = dd.from_pandas(df, npartitions=2)
    s1 = repr(ddf)
    assert repr(ddf) == s1
    ddf.b = ddf.b.astype('category')
    assert repr(ddf) != s1


def test_dataframe_format():
    df = pd.DataFrame({'A': [1, 2, 3, 4, 5, 6, 7, 8],
                       'B': list('ABCDEFGH'),
                       'C': pd.Categorical(list('AAABBBCC'))})
    ddf = dd.from_pandas(df, 3)
    exp = ("Dask DataFrame Structure:\n"
           "                   A       B                C\n"
           "npartitions=3                                \n"
           "0              int64  object  category[known]\n"
           "3                ...     ...              ...\n"
           "6                ...     ...              ...\n"
           "7                ...     ...              ...\n"
           "Dask Name: from_pandas, 3 tasks")
    assert repr(ddf) == exp
    assert str(ddf) == exp

    exp = ("                   A       B                C\n"
           "npartitions=3                                \n"
           "0              int64  object  category[known]\n"
           "3                ...     ...              ...\n"
           "6                ...     ...              ...\n"
           "7                ...     ...              ...")
    assert ddf.to_string() == exp

    exp_table = """<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>A</th>
      <th>B</th>
      <th>C</th>
    </tr>
    <tr>
      <th>npartitions=3</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>int64</td>
      <td>object</td>
      <td>category[known]</td>
    </tr>
    <tr>
      <th>3</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>6</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>7</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
  </tbody>
</table>"""

    exp = """<div><strong>Dask DataFrame Structure:</strong></div>
{exp_table}
<div>Dask Name: from_pandas, 3 tasks</div>""".format(exp_table=exp_table)
    assert ddf.to_html() == exp

    # table is boxed with div and has style
    exp = """<div><strong>Dask DataFrame Structure:</strong></div>
<div>
{style}{exp_table}
</div>
<div>Dask Name: from_pandas, 3 tasks</div>""".format(style=style, exp_table=exp_table)
    assert ddf._repr_html_() == exp


def test_dataframe_format_with_index():
    df = pd.DataFrame({'A': [1, 2, 3, 4, 5, 6, 7, 8],
                       'B': list('ABCDEFGH'),
                       'C': pd.Categorical(list('AAABBBCC'))},
                      index=list('ABCDEFGH'))
    ddf = dd.from_pandas(df, 3)
    exp = ("Dask DataFrame Structure:\n"
           "                   A       B                C\n"
           "npartitions=3                                \n"
           "A              int64  object  category[known]\n"
           "D                ...     ...              ...\n"
           "G                ...     ...              ...\n"
           "H                ...     ...              ...\n"
           "Dask Name: from_pandas, 3 tasks")
    assert repr(ddf) == exp
    assert str(ddf) == exp

    exp_table = """<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>A</th>
      <th>B</th>
      <th>C</th>
    </tr>
    <tr>
      <th>npartitions=3</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>A</th>
      <td>int64</td>
      <td>object</td>
      <td>category[known]</td>
    </tr>
    <tr>
      <th>D</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>G</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>H</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
  </tbody>
</table>"""

    exp = """<div><strong>Dask DataFrame Structure:</strong></div>
{exp_table}
<div>Dask Name: from_pandas, 3 tasks</div>""".format(exp_table=exp_table)
    assert ddf.to_html() == exp

    # table is boxed with div and has style
    exp = """<div><strong>Dask DataFrame Structure:</strong></div>
<div>
{style}{exp_table}
</div>
<div>Dask Name: from_pandas, 3 tasks</div>""".format(style=style, exp_table=exp_table)
    assert ddf._repr_html_() == exp


def test_dataframe_format_unknown_divisions():
    df = pd.DataFrame({'A': [1, 2, 3, 4, 5, 6, 7, 8],
                       'B': list('ABCDEFGH'),
                       'C': pd.Categorical(list('AAABBBCC'))})
    ddf = dd.from_pandas(df, 3)
    ddf = ddf.clear_divisions()
    assert not ddf.known_divisions

    exp = ("Dask DataFrame Structure:\n"
           "                   A       B                C\n"
           "npartitions=3                                \n"
           "               int64  object  category[known]\n"
           "                 ...     ...              ...\n"
           "                 ...     ...              ...\n"
           "                 ...     ...              ...\n"
           "Dask Name: from_pandas, 3 tasks")
    assert repr(ddf) == exp
    assert str(ddf) == exp

    exp = ("                   A       B                C\n"
           "npartitions=3                                \n"
           "               int64  object  category[known]\n"
           "                 ...     ...              ...\n"
           "                 ...     ...              ...\n"
           "                 ...     ...              ...")
    assert ddf.to_string() == exp

    exp_table = """<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>A</th>
      <th>B</th>
      <th>C</th>
    </tr>
    <tr>
      <th>npartitions=3</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th></th>
      <td>int64</td>
      <td>object</td>
      <td>category[known]</td>
    </tr>
    <tr>
      <th></th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th></th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th></th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
  </tbody>
</table>"""

    exp = """<div><strong>Dask DataFrame Structure:</strong></div>
{exp_table}
<div>Dask Name: from_pandas, 3 tasks</div>""".format(exp_table=exp_table)
    assert ddf.to_html() == exp

    # table is boxed with div and has style
    exp = """<div><strong>Dask DataFrame Structure:</strong></div>
<div>
{style}{exp_table}
</div>
<div>Dask Name: from_pandas, 3 tasks</div>""".format(style=style, exp_table=exp_table)
    assert ddf._repr_html_() == exp


def test_dataframe_format_long():
    df = pd.DataFrame({'A': [1, 2, 3, 4, 5, 6, 7, 8] * 10,
                       'B': list('ABCDEFGH') * 10,
                       'C': pd.Categorical(list('AAABBBCC') * 10)})
    ddf = dd.from_pandas(df, 10)
    exp = ('Dask DataFrame Structure:\n'
           '                    A       B                C\n'
           'npartitions=10                                \n'
           '0               int64  object  category[known]\n'
           '8                 ...     ...              ...\n'
           '...               ...     ...              ...\n'
           '72                ...     ...              ...\n'
           '79                ...     ...              ...\n'
           'Dask Name: from_pandas, 10 tasks')
    assert repr(ddf) == exp
    assert str(ddf) == exp

    exp = ("                    A       B                C\n"
           "npartitions=10                                \n"
           "0               int64  object  category[known]\n"
           "8                 ...     ...              ...\n"
           "...               ...     ...              ...\n"
           "72                ...     ...              ...\n"
           "79                ...     ...              ...")
    assert ddf.to_string() == exp

    exp_table = """<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>A</th>
      <th>B</th>
      <th>C</th>
    </tr>
    <tr>
      <th>npartitions=10</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>int64</td>
      <td>object</td>
      <td>category[known]</td>
    </tr>
    <tr>
      <th>8</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>72</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>79</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
  </tbody>
</table>"""

    exp = """<div><strong>Dask DataFrame Structure:</strong></div>
{exp_table}
<div>Dask Name: from_pandas, 10 tasks</div>""".format(exp_table=exp_table)
    assert ddf.to_html() == exp

    # table is boxed with div
    exp = u"""<div><strong>Dask DataFrame Structure:</strong></div>
<div>
{style}{exp_table}
</div>
<div>Dask Name: from_pandas, 10 tasks</div>""".format(style=style, exp_table=exp_table)
    assert ddf._repr_html_() == exp


def test_series_format():
    s = pd.Series([1, 2, 3, 4, 5, 6, 7, 8],
                  index=list('ABCDEFGH'))
    ds = dd.from_pandas(s, 3)
    exp = """Dask Series Structure:
npartitions=3
A    int64
D      ...
G      ...
H      ...
dtype: int64
Dask Name: from_pandas, 3 tasks"""
    assert repr(ds) == exp
    assert str(ds) == exp

    exp = """npartitions=3
A    int64
D      ...
G      ...
H      ..."""
    assert ds.to_string() == exp

    s = pd.Series([1, 2, 3, 4, 5, 6, 7, 8],
                  index=list('ABCDEFGH'), name='XXX')
    ds = dd.from_pandas(s, 3)
    exp = """Dask Series Structure:
npartitions=3
A    int64
D      ...
G      ...
H      ...
Name: XXX, dtype: int64
Dask Name: from_pandas, 3 tasks"""
    assert repr(ds) == exp
    assert str(ds) == exp


def test_series_format_long():
    s = pd.Series([1, 2, 3, 4, 5, 6, 7, 8, 9, 10] * 10,
                  index=list('ABCDEFGHIJ') * 10)
    ds = dd.from_pandas(s, 10)
    exp = ("Dask Series Structure:\nnpartitions=10\nA    int64\nB      ...\n"
           "     ...  \nJ      ...\nJ      ...\ndtype: int64\n"
           "Dask Name: from_pandas, 10 tasks")
    assert repr(ds) == exp
    assert str(ds) == exp

    exp = "npartitions=10\nA    int64\nB      ...\n     ...  \nJ      ...\nJ      ..."
    assert ds.to_string() == exp


def test_index_format():
    s = pd.Series([1, 2, 3, 4, 5, 6, 7, 8],
                  index=list('ABCDEFGH'))
    ds = dd.from_pandas(s, 3)
    exp = """Dask Index Structure:
npartitions=3
A    object
D       ...
G       ...
H       ...
dtype: object
Dask Name: from_pandas, 6 tasks"""
    assert repr(ds.index) == exp
    assert str(ds.index) == exp

    s = pd.Series([1, 2, 3, 4, 5, 6, 7, 8],
                  index=pd.CategoricalIndex([1, 2, 3, 4, 5, 6, 7, 8], name='YYY'))
    ds = dd.from_pandas(s, 3)
    exp = dedent("""\
    Dask Index Structure:
    npartitions=3
    1    category[known]
    4                ...
    7                ...
    8                ...
    Name: YYY, dtype: category
    Dask Name: from_pandas, 6 tasks""")
    assert repr(ds.index) == exp
    assert str(ds.index) == exp


def test_categorical_format():
    s = pd.Series(['a', 'b', 'c']).astype('category')
    known = dd.from_pandas(s, npartitions=1)
    unknown = known.cat.as_unknown()
    exp = ("Dask Series Structure:\n"
           "npartitions=1\n"
           "0    category[known]\n"
           "2                ...\n"
           "dtype: category\n"
           "Dask Name: from_pandas, 1 tasks")
    assert repr(known) == exp
    exp = ("Dask Series Structure:\n"
           "npartitions=1\n"
           "0    category[unknown]\n"
           "2                  ...\n"
           "dtype: category\n"
           "Dask Name: from_pandas, 1 tasks")
    assert repr(unknown) == exp
