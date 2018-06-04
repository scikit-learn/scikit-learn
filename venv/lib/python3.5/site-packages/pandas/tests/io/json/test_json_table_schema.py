"""Tests for Table Schema integration."""
import json
from collections import OrderedDict

import numpy as np
import pandas as pd
import pytest

from pandas import DataFrame
from pandas.core.dtypes.dtypes import (
    PeriodDtype, CategoricalDtype, DatetimeTZDtype)
from pandas.io.json.table_schema import (
    as_json_table_type,
    build_table_schema,
    convert_pandas_type_to_json_field,
    convert_json_field_to_pandas_type,
    set_default_names)
import pandas.util.testing as tm


class TestBuildSchema(object):

    def setup_method(self, method):
        self.df = DataFrame(
            {'A': [1, 2, 3, 4],
             'B': ['a', 'b', 'c', 'c'],
             'C': pd.date_range('2016-01-01', freq='d', periods=4),
             'D': pd.timedelta_range('1H', periods=4, freq='T'),
             },
            index=pd.Index(range(4), name='idx'))

    def test_build_table_schema(self):
        result = build_table_schema(self.df, version=False)
        expected = {
            'fields': [{'name': 'idx', 'type': 'integer'},
                       {'name': 'A', 'type': 'integer'},
                       {'name': 'B', 'type': 'string'},
                       {'name': 'C', 'type': 'datetime'},
                       {'name': 'D', 'type': 'duration'},
                       ],
            'primaryKey': ['idx']
        }
        assert result == expected
        result = build_table_schema(self.df)
        assert "pandas_version" in result

    def test_series(self):
        s = pd.Series([1, 2, 3], name='foo')
        result = build_table_schema(s, version=False)
        expected = {'fields': [{'name': 'index', 'type': 'integer'},
                               {'name': 'foo', 'type': 'integer'}],
                    'primaryKey': ['index']}
        assert result == expected
        result = build_table_schema(s)
        assert 'pandas_version' in result

    def test_series_unnamed(self):
        result = build_table_schema(pd.Series([1, 2, 3]), version=False)
        expected = {'fields': [{'name': 'index', 'type': 'integer'},
                               {'name': 'values', 'type': 'integer'}],
                    'primaryKey': ['index']}
        assert result == expected

    def test_multiindex(self):
        df = self.df.copy()
        idx = pd.MultiIndex.from_product([('a', 'b'), (1, 2)])
        df.index = idx

        result = build_table_schema(df, version=False)
        expected = {
            'fields': [{'name': 'level_0', 'type': 'string'},
                       {'name': 'level_1', 'type': 'integer'},
                       {'name': 'A', 'type': 'integer'},
                       {'name': 'B', 'type': 'string'},
                       {'name': 'C', 'type': 'datetime'},
                       {'name': 'D', 'type': 'duration'},
                       ],
            'primaryKey': ['level_0', 'level_1']
        }
        assert result == expected

        df.index.names = ['idx0', None]
        expected['fields'][0]['name'] = 'idx0'
        expected['primaryKey'] = ['idx0', 'level_1']
        result = build_table_schema(df, version=False)
        assert result == expected


class TestTableSchemaType(object):

    @pytest.mark.parametrize('int_type', [
        np.int, np.int16, np.int32, np.int64])
    def test_as_json_table_type_int_data(self, int_type):
        int_data = [1, 2, 3]
        assert as_json_table_type(np.array(
            int_data, dtype=int_type)) == 'integer'

    @pytest.mark.parametrize('float_type', [
        np.float, np.float16, np.float32, np.float64])
    def test_as_json_table_type_float_data(self, float_type):
        float_data = [1., 2., 3.]
        assert as_json_table_type(np.array(
            float_data, dtype=float_type)) == 'number'

    @pytest.mark.parametrize('bool_type', [bool, np.bool])
    def test_as_json_table_type_bool_data(self, bool_type):
        bool_data = [True, False]
        assert as_json_table_type(np.array(
            bool_data, dtype=bool_type)) == 'boolean'

    @pytest.mark.parametrize('date_data', [
        pd.to_datetime(['2016']),
        pd.to_datetime(['2016'], utc=True),
        pd.Series(pd.to_datetime(['2016'])),
        pd.Series(pd.to_datetime(['2016'], utc=True)),
        pd.period_range('2016', freq='A', periods=3)
    ])
    def test_as_json_table_type_date_data(self, date_data):
        assert as_json_table_type(date_data) == 'datetime'

    @pytest.mark.parametrize('str_data', [
        pd.Series(['a', 'b']), pd.Index(['a', 'b'])])
    def test_as_json_table_type_string_data(self, str_data):
        assert as_json_table_type(str_data) == 'string'

    @pytest.mark.parametrize('cat_data', [
        pd.Categorical(['a']),
        pd.Categorical([1]),
        pd.Series(pd.Categorical([1])),
        pd.CategoricalIndex([1]),
        pd.Categorical([1])])
    def test_as_json_table_type_categorical_data(self, cat_data):
        assert as_json_table_type(cat_data) == 'any'

    # ------
    # dtypes
    # ------
    @pytest.mark.parametrize('int_dtype', [
        np.int, np.int16, np.int32, np.int64])
    def test_as_json_table_type_int_dtypes(self, int_dtype):
        assert as_json_table_type(int_dtype) == 'integer'

    @pytest.mark.parametrize('float_dtype', [
        np.float, np.float16, np.float32, np.float64])
    def test_as_json_table_type_float_dtypes(self, float_dtype):
        assert as_json_table_type(float_dtype) == 'number'

    @pytest.mark.parametrize('bool_dtype', [bool, np.bool])
    def test_as_json_table_type_bool_dtypes(self, bool_dtype):
        assert as_json_table_type(bool_dtype) == 'boolean'

    @pytest.mark.parametrize('date_dtype', [
        np.datetime64, np.dtype("<M8[ns]"), PeriodDtype(),
        DatetimeTZDtype('ns', 'US/Central')])
    def test_as_json_table_type_date_dtypes(self, date_dtype):
        # TODO: datedate.date? datetime.time?
        assert as_json_table_type(date_dtype) == 'datetime'

    @pytest.mark.parametrize('td_dtype', [
        np.timedelta64, np.dtype("<m8[ns]")])
    def test_as_json_table_type_timedelta_dtypes(self, td_dtype):
        assert as_json_table_type(td_dtype) == 'duration'

    @pytest.mark.parametrize('str_dtype', [object])  # TODO
    def test_as_json_table_type_string_dtypes(self, str_dtype):
        assert as_json_table_type(str_dtype) == 'string'

    def test_as_json_table_type_categorical_dtypes(self):
        # TODO: I think before is_categorical_dtype(Categorical)
        # returned True, but now it's False. Figure out why or
        # if it matters
        assert as_json_table_type(pd.Categorical(['a'])) == 'any'
        assert as_json_table_type(CategoricalDtype()) == 'any'


class TestTableOrient(object):

    def setup_method(self, method):
        self.df = DataFrame(
            {'A': [1, 2, 3, 4],
             'B': ['a', 'b', 'c', 'c'],
             'C': pd.date_range('2016-01-01', freq='d', periods=4),
             'D': pd.timedelta_range('1H', periods=4, freq='T'),
             'E': pd.Series(pd.Categorical(['a', 'b', 'c', 'c'])),
             'F': pd.Series(pd.Categorical(['a', 'b', 'c', 'c'],
                                           ordered=True)),
             'G': [1., 2., 3, 4.],
             'H': pd.date_range('2016-01-01', freq='d', periods=4,
                                tz='US/Central'),
             },
            index=pd.Index(range(4), name='idx'))

    def test_build_series(self):
        s = pd.Series([1, 2], name='a')
        s.index.name = 'id'
        result = s.to_json(orient='table', date_format='iso')
        result = json.loads(result, object_pairs_hook=OrderedDict)

        assert "pandas_version" in result['schema']
        result['schema'].pop('pandas_version')

        fields = [{'name': 'id', 'type': 'integer'},
                  {'name': 'a', 'type': 'integer'}]

        schema = {
            'fields': fields,
            'primaryKey': ['id'],
        }

        expected = OrderedDict([
            ('schema', schema),
            ('data', [OrderedDict([('id', 0), ('a', 1)]),
                      OrderedDict([('id', 1), ('a', 2)])])])
        assert result == expected

    def test_to_json(self):
        df = self.df.copy()
        df.index.name = 'idx'
        result = df.to_json(orient='table', date_format='iso')
        result = json.loads(result, object_pairs_hook=OrderedDict)

        assert "pandas_version" in result['schema']
        result['schema'].pop('pandas_version')

        fields = [
            {'name': 'idx', 'type': 'integer'},
            {'name': 'A', 'type': 'integer'},
            {'name': 'B', 'type': 'string'},
            {'name': 'C', 'type': 'datetime'},
            {'name': 'D', 'type': 'duration'},
            {'constraints': {'enum': ['a', 'b', 'c']},
             'name': 'E',
             'ordered': False,
             'type': 'any'},
            {'constraints': {'enum': ['a', 'b', 'c']},
             'name': 'F',
             'ordered': True,
             'type': 'any'},
            {'name': 'G', 'type': 'number'},
            {'name': 'H', 'type': 'datetime', 'tz': 'US/Central'}
        ]

        schema = {
            'fields': fields,
            'primaryKey': ['idx'],
        }
        data = [
            OrderedDict([('idx', 0), ('A', 1), ('B', 'a'),
                         ('C', '2016-01-01T00:00:00.000Z'),
                         ('D', 'P0DT1H0M0S'),
                         ('E', 'a'), ('F', 'a'), ('G', 1.),
                         ('H', '2016-01-01T06:00:00.000Z')
                         ]),
            OrderedDict([('idx', 1), ('A', 2), ('B', 'b'),
                         ('C', '2016-01-02T00:00:00.000Z'),
                         ('D', 'P0DT1H1M0S'),
                         ('E', 'b'), ('F', 'b'), ('G', 2.),
                         ('H', '2016-01-02T06:00:00.000Z')
                         ]),
            OrderedDict([('idx', 2), ('A', 3), ('B', 'c'),
                         ('C', '2016-01-03T00:00:00.000Z'),
                         ('D', 'P0DT1H2M0S'),
                         ('E', 'c'), ('F', 'c'), ('G', 3.),
                         ('H', '2016-01-03T06:00:00.000Z')
                         ]),
            OrderedDict([('idx', 3), ('A', 4), ('B', 'c'),
                         ('C', '2016-01-04T00:00:00.000Z'),
                         ('D', 'P0DT1H3M0S'),
                         ('E', 'c'), ('F', 'c'), ('G', 4.),
                         ('H', '2016-01-04T06:00:00.000Z')
                         ]),
        ]
        expected = OrderedDict([('schema', schema), ('data', data)])
        assert result == expected

    def test_to_json_float_index(self):
        data = pd.Series(1, index=[1., 2.])
        result = data.to_json(orient='table', date_format='iso')
        result = json.loads(result, object_pairs_hook=OrderedDict)
        result['schema'].pop('pandas_version')

        expected = (
            OrderedDict([('schema', {
                'fields': [{'name': 'index', 'type': 'number'},
                           {'name': 'values', 'type': 'integer'}],
                'primaryKey': ['index']
            }),
                ('data', [OrderedDict([('index', 1.0), ('values', 1)]),
                          OrderedDict([('index', 2.0), ('values', 1)])])])
        )
        assert result == expected

    def test_to_json_period_index(self):
        idx = pd.period_range('2016', freq='Q-JAN', periods=2)
        data = pd.Series(1, idx)
        result = data.to_json(orient='table', date_format='iso')
        result = json.loads(result, object_pairs_hook=OrderedDict)
        result['schema'].pop('pandas_version')

        fields = [{'freq': 'Q-JAN', 'name': 'index', 'type': 'datetime'},
                  {'name': 'values', 'type': 'integer'}]

        schema = {'fields': fields, 'primaryKey': ['index']}
        data = [OrderedDict([('index', '2015-11-01T00:00:00.000Z'),
                             ('values', 1)]),
                OrderedDict([('index', '2016-02-01T00:00:00.000Z'),
                             ('values', 1)])]
        expected = OrderedDict([('schema', schema), ('data', data)])
        assert result == expected

    def test_to_json_categorical_index(self):
        data = pd.Series(1, pd.CategoricalIndex(['a', 'b']))
        result = data.to_json(orient='table', date_format='iso')
        result = json.loads(result, object_pairs_hook=OrderedDict)
        result['schema'].pop('pandas_version')

        expected = (
            OrderedDict([('schema',
                          {'fields': [{'name': 'index', 'type': 'any',
                                       'constraints': {'enum': ['a', 'b']},
                                       'ordered': False},
                                      {'name': 'values', 'type': 'integer'}],
                           'primaryKey': ['index']}),
                         ('data', [
                             OrderedDict([('index', 'a'),
                                          ('values', 1)]),
                             OrderedDict([('index', 'b'), ('values', 1)])])])
        )
        assert result == expected

    def test_date_format_raises(self):
        with pytest.raises(ValueError):
            self.df.to_json(orient='table', date_format='epoch')

        # others work
        self.df.to_json(orient='table', date_format='iso')
        self.df.to_json(orient='table')

    @pytest.mark.parametrize('kind', [pd.Series, pd.Index])
    def test_convert_pandas_type_to_json_field_int(self, kind):
        data = [1, 2, 3]
        result = convert_pandas_type_to_json_field(kind(data, name='name'))
        expected = {"name": "name", "type": "integer"}
        assert result == expected

    @pytest.mark.parametrize('kind', [pd.Series, pd.Index])
    def test_convert_pandas_type_to_json_field_float(self, kind):
        data = [1., 2., 3.]
        result = convert_pandas_type_to_json_field(kind(data, name='name'))
        expected = {"name": "name", "type": "number"}
        assert result == expected

    @pytest.mark.parametrize('dt_args,extra_exp', [
        ({}, {}), ({'utc': True}, {'tz': 'UTC'})])
    @pytest.mark.parametrize('wrapper', [None, pd.Series])
    def test_convert_pandas_type_to_json_field_datetime(self, dt_args,
                                                        extra_exp, wrapper):
        data = [1., 2., 3.]
        data = pd.to_datetime(data, **dt_args)
        if wrapper is pd.Series:
            data = pd.Series(data, name='values')
        result = convert_pandas_type_to_json_field(data)
        expected = {"name": "values", "type": 'datetime'}
        expected.update(extra_exp)
        assert result == expected

    def test_convert_pandas_type_to_json_period_range(self):
        arr = pd.period_range('2016', freq='A-DEC', periods=4)
        result = convert_pandas_type_to_json_field(arr)
        expected = {"name": "values", "type": 'datetime', "freq": "A-DEC"}
        assert result == expected

    @pytest.mark.parametrize('kind', [pd.Categorical, pd.CategoricalIndex])
    @pytest.mark.parametrize('ordered', [True, False])
    def test_convert_pandas_type_to_json_field_categorical(self, kind,
                                                           ordered):
        data = ['a', 'b', 'c']
        if kind is pd.Categorical:
            arr = pd.Series(kind(data, ordered=ordered), name='cats')
        elif kind is pd.CategoricalIndex:
            arr = kind(data, ordered=ordered, name='cats')

        result = convert_pandas_type_to_json_field(arr)
        expected = {"name": "cats", "type": "any",
                    "constraints": {"enum": data},
                    "ordered": ordered}
        assert result == expected

    @pytest.mark.parametrize("inp,exp", [
        ({'type': 'integer'}, 'int64'),
        ({'type': 'number'}, 'float64'),
        ({'type': 'boolean'}, 'bool'),
        ({'type': 'duration'}, 'timedelta64'),
        ({'type': 'datetime'}, 'datetime64[ns]'),
        ({'type': 'datetime', 'tz': 'US/Hawaii'}, 'datetime64[ns, US/Hawaii]'),
        ({'type': 'any'}, 'object'),
        ({'type': 'any', 'constraints': {'enum': ['a', 'b', 'c']},
          'ordered': False}, CategoricalDtype(categories=['a', 'b', 'c'],
                                              ordered=False)),
        ({'type': 'any', 'constraints': {'enum': ['a', 'b', 'c']},
          'ordered': True}, CategoricalDtype(categories=['a', 'b', 'c'],
                                             ordered=True)),
        ({'type': 'string'}, 'object')])
    def test_convert_json_field_to_pandas_type(self, inp, exp):
        field = {'name': 'foo'}
        field.update(inp)
        assert convert_json_field_to_pandas_type(field) == exp

    @pytest.mark.parametrize("inp", ["geopoint", "geojson", "fake_type"])
    def test_convert_json_field_to_pandas_type_raises(self, inp):
        field = {'type': inp}
        with tm.assert_raises_regex(ValueError, "Unsupported or invalid field "
                                    "type: {}".format(inp)):
            convert_json_field_to_pandas_type(field)

    def test_categorical(self):
        s = pd.Series(pd.Categorical(['a', 'b', 'a']))
        s.index.name = 'idx'
        result = s.to_json(orient='table', date_format='iso')
        result = json.loads(result, object_pairs_hook=OrderedDict)
        result['schema'].pop('pandas_version')

        fields = [{'name': 'idx', 'type': 'integer'},
                  {'constraints': {'enum': ['a', 'b']},
                   'name': 'values',
                   'ordered': False,
                   'type': 'any'}]

        expected = OrderedDict([
            ('schema', {'fields': fields,
                        'primaryKey': ['idx']}),
            ('data', [OrderedDict([('idx', 0), ('values', 'a')]),
                      OrderedDict([('idx', 1), ('values', 'b')]),
                      OrderedDict([('idx', 2), ('values', 'a')])])])
        assert result == expected

    @pytest.mark.parametrize('idx,nm,prop', [
        (pd.Index([1]), 'index', 'name'),
        (pd.Index([1], name='myname'), 'myname', 'name'),
        (pd.MultiIndex.from_product([('a', 'b'), ('c', 'd')]),
         ['level_0', 'level_1'], 'names'),
        (pd.MultiIndex.from_product([('a', 'b'), ('c', 'd')],
                                    names=['n1', 'n2']),
         ['n1', 'n2'], 'names'),
        (pd.MultiIndex.from_product([('a', 'b'), ('c', 'd')],
                                    names=['n1', None]),
         ['n1', 'level_1'], 'names')
    ])
    def test_set_names_unset(self, idx, nm, prop):
        data = pd.Series(1, idx)
        result = set_default_names(data)
        assert getattr(result.index, prop) == nm

    @pytest.mark.parametrize("idx", [
        pd.Index([], name='index'),
        pd.MultiIndex.from_arrays([['foo'], ['bar']],
                                  names=('level_0', 'level_1')),
        pd.MultiIndex.from_arrays([['foo'], ['bar']],
                                  names=('foo', 'level_1'))
    ])
    def test_warns_non_roundtrippable_names(self, idx):
        # GH 19130
        df = pd.DataFrame([[]], index=idx)
        df.index.name = 'index'
        with tm.assert_produces_warning():
            set_default_names(df)

    def test_timestamp_in_columns(self):
        df = pd.DataFrame([[1, 2]], columns=[pd.Timestamp('2016'),
                                             pd.Timedelta(10, unit='s')])
        result = df.to_json(orient="table")
        js = json.loads(result)
        assert js['schema']['fields'][1]['name'] == 1451606400000
        assert js['schema']['fields'][2]['name'] == 10000

    @pytest.mark.parametrize('case', [
        pd.Series([1], index=pd.Index([1], name='a'), name='a'),
        pd.DataFrame({"A": [1]}, index=pd.Index([1], name="A")),
        pd.DataFrame({"A": [1]}, index=pd.MultiIndex.from_arrays([
            ['a'], [1]], names=["A", "a"]))
    ])
    def test_overlapping_names(self, case):
        with tm.assert_raises_regex(ValueError, 'Overlapping'):
            case.to_json(orient='table')

    def test_mi_falsey_name(self):
        # GH 16203
        df = pd.DataFrame(np.random.randn(4, 4),
                          index=pd.MultiIndex.from_product([('A', 'B'),
                                                            ('a', 'b')]))
        result = [x['name'] for x in build_table_schema(df)['fields']]
        assert result == ['level_0', 'level_1', 0, 1, 2, 3]


class TestTableOrientReader(object):

    @pytest.mark.parametrize("index_nm", [
        None, "idx", pytest.param("index", marks=pytest.mark.xfail),
        'level_0'])
    @pytest.mark.parametrize("vals", [
        {'ints': [1, 2, 3, 4]},
        {'objects': ['a', 'b', 'c', 'd']},
        {'date_ranges': pd.date_range('2016-01-01', freq='d', periods=4)},
        {'categoricals': pd.Series(pd.Categorical(['a', 'b', 'c', 'c']))},
        {'ordered_cats': pd.Series(pd.Categorical(['a', 'b', 'c', 'c'],
                                                  ordered=True))},
        pytest.param({'floats': [1., 2., 3., 4.]}, marks=pytest.mark.xfail),
        {'floats': [1.1, 2.2, 3.3, 4.4]},
        {'bools': [True, False, False, True]}])
    def test_read_json_table_orient(self, index_nm, vals, recwarn):
        df = DataFrame(vals, index=pd.Index(range(4), name=index_nm))
        out = df.to_json(orient="table")
        result = pd.read_json(out, orient="table")
        tm.assert_frame_equal(df, result)

    @pytest.mark.parametrize("index_nm", [
        None, "idx", "index"])
    @pytest.mark.parametrize("vals", [
        {'timedeltas': pd.timedelta_range('1H', periods=4, freq='T')},
        {'timezones': pd.date_range('2016-01-01', freq='d', periods=4,
                                    tz='US/Central')}])
    def test_read_json_table_orient_raises(self, index_nm, vals, recwarn):
        df = DataFrame(vals, index=pd.Index(range(4), name=index_nm))
        out = df.to_json(orient="table")
        with tm.assert_raises_regex(NotImplementedError, 'can not yet read '):
            pd.read_json(out, orient="table")

    def test_comprehensive(self):
        df = DataFrame(
            {'A': [1, 2, 3, 4],
             'B': ['a', 'b', 'c', 'c'],
             'C': pd.date_range('2016-01-01', freq='d', periods=4),
             # 'D': pd.timedelta_range('1H', periods=4, freq='T'),
             'E': pd.Series(pd.Categorical(['a', 'b', 'c', 'c'])),
             'F': pd.Series(pd.Categorical(['a', 'b', 'c', 'c'],
                                           ordered=True)),
             'G': [1.1, 2.2, 3.3, 4.4],
             # 'H': pd.date_range('2016-01-01', freq='d', periods=4,
             #                   tz='US/Central'),
             'I': [True, False, False, True],
             },
            index=pd.Index(range(4), name='idx'))

        out = df.to_json(orient="table")
        result = pd.read_json(out, orient="table")
        tm.assert_frame_equal(df, result)

    @pytest.mark.parametrize("index_names", [
        [None, None], ['foo', 'bar'], ['foo', None], [None, 'foo'],
        ['index', 'foo']])
    def test_multiindex(self, index_names):
        # GH 18912
        df = pd.DataFrame(
            [["Arr", "alpha", [1, 2, 3, 4]],
             ["Bee", "Beta", [10, 20, 30, 40]]],
            index=[["A", "B"], ["Null", "Eins"]],
            columns=["Aussprache", "Griechisch", "Args"]
        )
        df.index.names = index_names
        out = df.to_json(orient="table")
        result = pd.read_json(out, orient="table")
        tm.assert_frame_equal(df, result)
