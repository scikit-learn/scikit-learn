from __future__ import annotations

from pandas._typing import ReadBuffer
from pandas.compat._optional import import_optional_dependency

from pandas.core.dtypes.inference import is_integer

from pandas.core.frame import DataFrame

from pandas.io.parsers.base_parser import ParserBase


class ArrowParserWrapper(ParserBase):
    """
    Wrapper for the pyarrow engine for read_csv()
    """

    def __init__(self, src: ReadBuffer[bytes], **kwds):
        super().__init__(kwds)
        self.kwds = kwds
        self.src = src

        self._parse_kwds()

    def _parse_kwds(self):
        """
        Validates keywords before passing to pyarrow.
        """
        encoding: str | None = self.kwds.get("encoding")
        self.encoding = "utf-8" if encoding is None else encoding

        self.usecols, self.usecols_dtype = self._validate_usecols_arg(
            self.kwds["usecols"]
        )
        na_values = self.kwds["na_values"]
        if isinstance(na_values, dict):
            raise ValueError(
                "The pyarrow engine doesn't support passing a dict for na_values"
            )
        self.na_values = list(self.kwds["na_values"])

    def _get_pyarrow_options(self):
        """
        Rename some arguments to pass to pyarrow
        """
        mapping = {
            "usecols": "include_columns",
            "na_values": "null_values",
            "escapechar": "escape_char",
            "skip_blank_lines": "ignore_empty_lines",
        }
        for pandas_name, pyarrow_name in mapping.items():
            if pandas_name in self.kwds and self.kwds.get(pandas_name) is not None:
                self.kwds[pyarrow_name] = self.kwds.pop(pandas_name)

        self.parse_options = {
            option_name: option_value
            for option_name, option_value in self.kwds.items()
            if option_value is not None
            and option_name
            in ("delimiter", "quote_char", "escape_char", "ignore_empty_lines")
        }
        self.convert_options = {
            option_name: option_value
            for option_name, option_value in self.kwds.items()
            if option_value is not None
            and option_name
            in ("include_columns", "null_values", "true_values", "false_values")
        }
        self.read_options = {
            "autogenerate_column_names": self.header is None,
            "skip_rows": self.header
            if self.header is not None
            else self.kwds["skiprows"],
        }

    def _finalize_output(self, frame: DataFrame) -> DataFrame:
        """
        Processes data read in based on kwargs.

        Parameters
        ----------
        frame: DataFrame
            The DataFrame to process.

        Returns
        -------
        DataFrame
            The processed DataFrame.
        """
        num_cols = len(frame.columns)
        multi_index_named = True
        if self.header is None:
            if self.names is None:
                if self.prefix is not None:
                    self.names = [f"{self.prefix}{i}" for i in range(num_cols)]
                elif self.header is None:
                    self.names = range(num_cols)
            if len(self.names) != num_cols:
                # usecols is passed through to pyarrow, we only handle index col here
                # The only way self.names is not the same length as number of cols is
                # if we have int index_col. We should just pad the names(they will get
                # removed anyways) to expected length then.
                self.names = list(range(num_cols - len(self.names))) + self.names
                multi_index_named = False
            frame.columns = self.names
        # we only need the frame not the names
        # error: Incompatible types in assignment (expression has type
        # "Union[List[Union[Union[str, int, float, bool], Union[Period, Timestamp,
        # Timedelta, Any]]], Index]", variable has type "Index")  [assignment]
        frame.columns, frame = self._do_date_conversions(  # type: ignore[assignment]
            frame.columns, frame
        )
        if self.index_col is not None:
            for i, item in enumerate(self.index_col):
                if is_integer(item):
                    self.index_col[i] = frame.columns[item]
                else:
                    # String case
                    if item not in frame.columns:
                        raise ValueError(f"Index {item} invalid")
            frame.set_index(self.index_col, drop=True, inplace=True)
            # Clear names if headerless and no name given
            if self.header is None and not multi_index_named:
                frame.index.names = [None] * len(frame.index.names)

        if self.kwds.get("dtype") is not None:
            try:
                frame = frame.astype(self.kwds.get("dtype"))
            except TypeError as e:
                # GH#44901 reraise to keep api consistent
                raise ValueError(e)
        return frame

    def read(self) -> DataFrame:
        """
        Reads the contents of a CSV file into a DataFrame and
        processes it according to the kwargs passed in the
        constructor.

        Returns
        -------
        DataFrame
            The DataFrame created from the CSV file.
        """
        pyarrow_csv = import_optional_dependency("pyarrow.csv")
        self._get_pyarrow_options()

        table = pyarrow_csv.read_csv(
            self.src,
            read_options=pyarrow_csv.ReadOptions(**self.read_options),
            parse_options=pyarrow_csv.ParseOptions(**self.parse_options),
            convert_options=pyarrow_csv.ConvertOptions(**self.convert_options),
        )

        frame = table.to_pandas()
        return self._finalize_output(frame)
