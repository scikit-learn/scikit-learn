import re

import pytest

import pandas._testing as tm

from pandas.io.excel import ExcelWriter

odf = pytest.importorskip("odf")

pytestmark = pytest.mark.parametrize("ext", [".ods"])


def test_write_append_mode_raises(ext):
    msg = "Append mode is not supported with odf!"

    with tm.ensure_clean(ext) as f:
        with pytest.raises(ValueError, match=msg):
            ExcelWriter(f, engine="odf", mode="a")


def test_kwargs(ext):
    # GH 42286
    # GH 43445
    # test for error: OpenDocumentSpreadsheet does not accept any arguments
    kwargs = {"kwarg": 1}
    with tm.ensure_clean(ext) as f:
        msg = re.escape("Use of **kwargs is deprecated")
        error = re.escape(
            "OpenDocumentSpreadsheet() got an unexpected keyword argument 'kwarg'"
        )
        with pytest.raises(
            TypeError,
            match=error,
        ):
            with tm.assert_produces_warning(FutureWarning, match=msg):
                with ExcelWriter(f, engine="odf", **kwargs) as _:
                    pass


@pytest.mark.parametrize("engine_kwargs", [None, {"kwarg": 1}])
def test_engine_kwargs(ext, engine_kwargs):
    # GH 42286
    # GH 43445
    # test for error: OpenDocumentSpreadsheet does not accept any arguments
    with tm.ensure_clean(ext) as f:
        if engine_kwargs is not None:
            error = re.escape(
                "OpenDocumentSpreadsheet() got an unexpected keyword argument 'kwarg'"
            )
            with pytest.raises(
                TypeError,
                match=error,
            ):
                ExcelWriter(f, engine="odf", engine_kwargs=engine_kwargs)
        else:
            with ExcelWriter(f, engine="odf", engine_kwargs=engine_kwargs) as _:
                pass
