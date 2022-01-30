import numpy as np
import pytest

from pandas import DataFrame
import pandas._testing as tm

from pandas.io.excel import ExcelWriter
from pandas.io.formats.excel import ExcelFormatter

pytest.importorskip("jinja2")
# jinja2 is currently required for Styler.__init__(). Technically Styler.to_excel
# could compute styles and render to excel without jinja2, since there is no
# 'template' file, but this needs the import error to delayed until render time.


def assert_equal_cell_styles(cell1, cell2):
    # TODO: should find a better way to check equality
    assert cell1.alignment.__dict__ == cell2.alignment.__dict__
    assert cell1.border.__dict__ == cell2.border.__dict__
    assert cell1.fill.__dict__ == cell2.fill.__dict__
    assert cell1.font.__dict__ == cell2.font.__dict__
    assert cell1.number_format == cell2.number_format
    assert cell1.protection.__dict__ == cell2.protection.__dict__


@pytest.mark.parametrize(
    "engine",
    ["xlsxwriter", "openpyxl"],
)
def test_styler_to_excel_unstyled(engine):
    # compare DataFrame.to_excel and Styler.to_excel when no styles applied
    pytest.importorskip(engine)
    df = DataFrame(np.random.randn(2, 2))
    with tm.ensure_clean(".xlsx") as path:
        with ExcelWriter(path, engine=engine) as writer:
            df.to_excel(writer, sheet_name="dataframe")
            df.style.to_excel(writer, sheet_name="unstyled")

        openpyxl = pytest.importorskip("openpyxl")  # test loading only with openpyxl
        wb = openpyxl.load_workbook(path)

        for col1, col2 in zip(wb["dataframe"].columns, wb["unstyled"].columns):
            assert len(col1) == len(col2)
            for cell1, cell2 in zip(col1, col2):
                assert cell1.value == cell2.value
                assert_equal_cell_styles(cell1, cell2)


shared_style_params = [
    (
        "background-color: #111222",
        ["fill", "fgColor", "rgb"],
        {"xlsxwriter": "FF111222", "openpyxl": "00111222"},
    ),
    (
        "color: #111222",
        ["font", "color", "value"],
        {"xlsxwriter": "FF111222", "openpyxl": "00111222"},
    ),
    ("font-family: Arial;", ["font", "name"], "arial"),
    ("font-weight: bold;", ["font", "b"], True),
    ("font-style: italic;", ["font", "i"], True),
    ("text-decoration: underline;", ["font", "u"], "single"),
    ("number-format: $??,???.00;", ["number_format"], "$??,???.00"),
    ("text-align: left;", ["alignment", "horizontal"], "left"),
    (
        "vertical-align: bottom;",
        ["alignment", "vertical"],
        {"xlsxwriter": None, "openpyxl": "bottom"},  # xlsxwriter Fails
    ),
]


@pytest.mark.parametrize(
    "engine",
    ["xlsxwriter", "openpyxl"],
)
@pytest.mark.parametrize("css, attrs, expected", shared_style_params)
def test_styler_to_excel_basic(engine, css, attrs, expected):
    pytest.importorskip(engine)
    df = DataFrame(np.random.randn(1, 1))
    styler = df.style.applymap(lambda x: css)

    with tm.ensure_clean(".xlsx") as path:
        with ExcelWriter(path, engine=engine) as writer:
            df.to_excel(writer, sheet_name="dataframe")
            styler.to_excel(writer, sheet_name="styled")

        openpyxl = pytest.importorskip("openpyxl")  # test loading only with openpyxl
        wb = openpyxl.load_workbook(path)

        # test unstyled data cell does not have expected styles
        # test styled cell has expected styles
        u_cell, s_cell = wb["dataframe"].cell(2, 2), wb["styled"].cell(2, 2)
        for attr in attrs:
            u_cell, s_cell = getattr(u_cell, attr), getattr(s_cell, attr)

        if isinstance(expected, dict):
            assert u_cell is None or u_cell != expected[engine]
            assert s_cell == expected[engine]
        else:
            assert u_cell is None or u_cell != expected
            assert s_cell == expected


@pytest.mark.parametrize(
    "engine",
    ["xlsxwriter", "openpyxl"],
)
@pytest.mark.parametrize("css, attrs, expected", shared_style_params)
def test_styler_to_excel_basic_indexes(engine, css, attrs, expected):
    pytest.importorskip(engine)
    df = DataFrame(np.random.randn(1, 1))

    styler = df.style
    styler.applymap_index(lambda x: css, axis=0)
    styler.applymap_index(lambda x: css, axis=1)

    null_styler = df.style
    null_styler.applymap(lambda x: "null: css;")
    null_styler.applymap_index(lambda x: "null: css;", axis=0)
    null_styler.applymap_index(lambda x: "null: css;", axis=1)

    with tm.ensure_clean(".xlsx") as path:
        with ExcelWriter(path, engine=engine) as writer:
            null_styler.to_excel(writer, sheet_name="null_styled")
            styler.to_excel(writer, sheet_name="styled")

        openpyxl = pytest.importorskip("openpyxl")  # test loading only with openpyxl
        wb = openpyxl.load_workbook(path)

        # test null styled index cells does not have expected styles
        # test styled cell has expected styles
        ui_cell, si_cell = wb["null_styled"].cell(2, 1), wb["styled"].cell(2, 1)
        uc_cell, sc_cell = wb["null_styled"].cell(1, 2), wb["styled"].cell(1, 2)
        for attr in attrs:
            ui_cell, si_cell = getattr(ui_cell, attr), getattr(si_cell, attr)
            uc_cell, sc_cell = getattr(uc_cell, attr), getattr(sc_cell, attr)

        if isinstance(expected, dict):
            assert ui_cell is None or ui_cell != expected[engine]
            assert si_cell == expected[engine]
            assert uc_cell is None or uc_cell != expected[engine]
            assert sc_cell == expected[engine]
        else:
            assert ui_cell is None or ui_cell != expected
            assert si_cell == expected
            assert uc_cell is None or uc_cell != expected
            assert sc_cell == expected


def test_styler_custom_converter():
    openpyxl = pytest.importorskip("openpyxl")

    def custom_converter(css):
        return {"font": {"color": {"rgb": "111222"}}}

    df = DataFrame(np.random.randn(1, 1))
    styler = df.style.applymap(lambda x: "color: #888999")
    with tm.ensure_clean(".xlsx") as path:
        with ExcelWriter(path, engine="openpyxl") as writer:
            ExcelFormatter(styler, style_converter=custom_converter).write(
                writer, sheet_name="custom"
            )

        wb = openpyxl.load_workbook(path)
        assert wb["custom"].cell(2, 2).font.color.value == "00111222"
