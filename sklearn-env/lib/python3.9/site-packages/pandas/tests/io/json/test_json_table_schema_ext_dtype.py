"""Tests for ExtensionDtype Table Schema integration."""

from collections import OrderedDict
import datetime as dt
import decimal
import json

import pytest

from pandas import (
    DataFrame,
    array,
)
from pandas.core.arrays.integer import Int64Dtype
from pandas.core.arrays.string_ import StringDtype
from pandas.core.series import Series
from pandas.tests.extension.date import (
    DateArray,
    DateDtype,
)
from pandas.tests.extension.decimal.array import (
    DecimalArray,
    DecimalDtype,
)

from pandas.io.json._table_schema import (
    as_json_table_type,
    build_table_schema,
)


class TestBuildSchema:
    def setup_method(self, method):
        self.da = DateArray([dt.date(2021, 10, 10)])
        self.dc = DecimalArray([decimal.Decimal(10)])
        self.sa = array(["pandas"], dtype="string")
        self.ia = array([10], dtype="Int64")
        self.df = DataFrame(
            {
                "A": self.da,
                "B": self.dc,
                "C": self.sa,
                "D": self.ia,
            }
        )

    def test_build_table_schema(self):
        result = build_table_schema(self.df, version=False)
        expected = {
            "fields": [
                {"name": "index", "type": "integer"},
                {"name": "A", "type": "any", "extDtype": "DateDtype"},
                {"name": "B", "type": "any", "extDtype": "decimal"},
                {"name": "C", "type": "any", "extDtype": "string"},
                {"name": "D", "type": "integer", "extDtype": "Int64"},
            ],
            "primaryKey": ["index"],
        }
        assert result == expected
        result = build_table_schema(self.df)
        assert "pandas_version" in result


class TestTableSchemaType:
    @pytest.mark.parametrize(
        "date_data",
        [
            DateArray([dt.date(2021, 10, 10)]),
            DateArray(dt.date(2021, 10, 10)),
            Series(DateArray(dt.date(2021, 10, 10))),
        ],
    )
    def test_as_json_table_type_ext_date_array_dtype(self, date_data):
        assert as_json_table_type(date_data.dtype) == "any"

    def test_as_json_table_type_ext_date_dtype(self):
        assert as_json_table_type(DateDtype()) == "any"

    @pytest.mark.parametrize(
        "decimal_data",
        [
            DecimalArray([decimal.Decimal(10)]),
            Series(DecimalArray([decimal.Decimal(10)])),
        ],
    )
    def test_as_json_table_type_ext_decimal_array_dtype(self, decimal_data):
        assert as_json_table_type(decimal_data.dtype) == "any"

    def test_as_json_table_type_ext_decimal_dtype(self):
        assert as_json_table_type(DecimalDtype()) == "any"

    @pytest.mark.parametrize(
        "string_data",
        [
            array(["pandas"], dtype="string"),
            Series(array(["pandas"], dtype="string")),
        ],
    )
    def test_as_json_table_type_ext_string_array_dtype(self, string_data):
        assert as_json_table_type(string_data.dtype) == "any"

    def test_as_json_table_type_ext_string_dtype(self):
        assert as_json_table_type(StringDtype()) == "any"

    @pytest.mark.parametrize(
        "integer_data",
        [
            array([10], dtype="Int64"),
            Series(array([10], dtype="Int64")),
        ],
    )
    def test_as_json_table_type_ext_integer_array_dtype(self, integer_data):
        assert as_json_table_type(integer_data.dtype) == "integer"

    def test_as_json_table_type_ext_integer_dtype(self):
        assert as_json_table_type(Int64Dtype()) == "integer"


class TestTableOrient:
    def setup_method(self, method):
        self.da = DateArray([dt.date(2021, 10, 10)])
        self.dc = DecimalArray([decimal.Decimal(10)])
        self.sa = array(["pandas"], dtype="string")
        self.ia = array([10], dtype="Int64")
        self.df = DataFrame(
            {
                "A": self.da,
                "B": self.dc,
                "C": self.sa,
                "D": self.ia,
            }
        )

    def test_build_date_series(self):

        s = Series(self.da, name="a")
        s.index.name = "id"
        result = s.to_json(orient="table", date_format="iso")
        result = json.loads(result, object_pairs_hook=OrderedDict)

        assert "pandas_version" in result["schema"]
        result["schema"].pop("pandas_version")

        fields = [
            {"name": "id", "type": "integer"},
            {"name": "a", "type": "any", "extDtype": "DateDtype"},
        ]

        schema = {"fields": fields, "primaryKey": ["id"]}

        expected = OrderedDict(
            [
                ("schema", schema),
                ("data", [OrderedDict([("id", 0), ("a", "2021-10-10T00:00:00.000Z")])]),
            ]
        )

        assert result == expected

    def test_build_decimal_series(self):

        s = Series(self.dc, name="a")
        s.index.name = "id"
        result = s.to_json(orient="table", date_format="iso")
        result = json.loads(result, object_pairs_hook=OrderedDict)

        assert "pandas_version" in result["schema"]
        result["schema"].pop("pandas_version")

        fields = [
            {"name": "id", "type": "integer"},
            {"name": "a", "type": "any", "extDtype": "decimal"},
        ]

        schema = {"fields": fields, "primaryKey": ["id"]}

        expected = OrderedDict(
            [
                ("schema", schema),
                ("data", [OrderedDict([("id", 0), ("a", 10.0)])]),
            ]
        )

        assert result == expected

    def test_build_string_series(self):
        s = Series(self.sa, name="a")
        s.index.name = "id"
        result = s.to_json(orient="table", date_format="iso")
        result = json.loads(result, object_pairs_hook=OrderedDict)

        assert "pandas_version" in result["schema"]
        result["schema"].pop("pandas_version")

        fields = [
            {"name": "id", "type": "integer"},
            {"name": "a", "type": "any", "extDtype": "string"},
        ]

        schema = {"fields": fields, "primaryKey": ["id"]}

        expected = OrderedDict(
            [
                ("schema", schema),
                ("data", [OrderedDict([("id", 0), ("a", "pandas")])]),
            ]
        )

        assert result == expected

    def test_build_int64_series(self):
        s = Series(self.ia, name="a")
        s.index.name = "id"
        result = s.to_json(orient="table", date_format="iso")
        result = json.loads(result, object_pairs_hook=OrderedDict)

        assert "pandas_version" in result["schema"]
        result["schema"].pop("pandas_version")

        fields = [
            {"name": "id", "type": "integer"},
            {"name": "a", "type": "integer", "extDtype": "Int64"},
        ]

        schema = {"fields": fields, "primaryKey": ["id"]}

        expected = OrderedDict(
            [
                ("schema", schema),
                ("data", [OrderedDict([("id", 0), ("a", 10)])]),
            ]
        )

        assert result == expected

    def test_to_json(self):

        df = self.df.copy()
        df.index.name = "idx"
        result = df.to_json(orient="table", date_format="iso")
        result = json.loads(result, object_pairs_hook=OrderedDict)

        assert "pandas_version" in result["schema"]
        result["schema"].pop("pandas_version")

        fields = [
            OrderedDict({"name": "idx", "type": "integer"}),
            OrderedDict({"name": "A", "type": "any", "extDtype": "DateDtype"}),
            OrderedDict({"name": "B", "type": "any", "extDtype": "decimal"}),
            OrderedDict({"name": "C", "type": "any", "extDtype": "string"}),
            OrderedDict({"name": "D", "type": "integer", "extDtype": "Int64"}),
        ]

        schema = OrderedDict({"fields": fields, "primaryKey": ["idx"]})
        data = [
            OrderedDict(
                [
                    ("idx", 0),
                    ("A", "2021-10-10T00:00:00.000Z"),
                    ("B", 10.0),
                    ("C", "pandas"),
                    ("D", 10),
                ]
            )
        ]
        expected = OrderedDict([("schema", schema), ("data", data)])

        assert result == expected
