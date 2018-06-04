# -*- coding: utf-8 -*-
"""
Internal module for formatting output data in csv, html,
and latex files. This module also applies to display formatting.
"""

from __future__ import print_function
# pylint: disable=W0141

from pandas.core.dtypes.missing import isna, notna
from pandas.core.dtypes.common import (
    is_categorical_dtype,
    is_float_dtype,
    is_period_arraylike,
    is_integer_dtype,
    is_interval_dtype,
    is_datetimetz,
    is_integer,
    is_float,
    is_scalar,
    is_numeric_dtype,
    is_datetime64_dtype,
    is_timedelta64_dtype,
    is_list_like)
from pandas.core.dtypes.generic import ABCSparseArray
from pandas.core.base import PandasObject
import pandas.core.common as com
from pandas.core.index import Index, MultiIndex, _ensure_index
from pandas import compat
from pandas.compat import (StringIO, lzip, map, zip, u)

from pandas.io.formats.terminal import get_terminal_size
from pandas.core.config import get_option, set_option
from pandas.io.common import (_expand_user, _stringify_path)
from pandas.io.formats.printing import adjoin, justify, pprint_thing
from pandas._libs import lib

from pandas._libs.tslib import (iNaT, Timestamp, Timedelta,
                                format_array_from_datetime)
from pandas.core.indexes.datetimes import DatetimeIndex
from pandas.core.indexes.period import PeriodIndex
import pandas as pd
import numpy as np

from functools import partial

common_docstring = """
    Parameters
    ----------
    buf : StringIO-like, optional
        buffer to write to
    columns : sequence, optional
        the subset of columns to write; default None writes all columns
    col_space : int, optional
        the minimum width of each column
    header : bool, optional
        %(header)s
    index : bool, optional
        whether to print index (row) labels, default True
    na_rep : string, optional
        string representation of NAN to use, default 'NaN'
    formatters : list or dict of one-parameter functions, optional
        formatter functions to apply to columns' elements by position or name,
        default None. The result of each function must be a unicode string.
        List must be of length equal to the number of columns.
    float_format : one-parameter function, optional
        formatter function to apply to columns' elements if they are floats,
        default None. The result of this function must be a unicode string.
    sparsify : bool, optional
        Set to False for a DataFrame with a hierarchical index to print every
        multiindex key at each row, default True
    index_names : bool, optional
        Prints the names of the indexes, default True
    line_width : int, optional
        Width to wrap a line in characters, default no wrap
    table_id : str, optional
        id for the <table> element create by to_html

        .. versionadded:: 0.23.0"""

_VALID_JUSTIFY_PARAMETERS = ("left", "right", "center", "justify",
                             "justify-all", "start", "end", "inherit",
                             "match-parent", "initial", "unset")

justify_docstring = """
    justify : str, default None
        How to justify the column labels. If None uses the option from
        the print configuration (controlled by set_option), 'right' out
        of the box. Valid values are

        * left
        * right
        * center
        * justify
        * justify-all
        * start
        * end
        * inherit
        * match-parent
        * initial
        * unset
"""

return_docstring = """

    Returns
    -------
    formatted : string (or unicode, depending on data and options)"""

docstring_to_string = common_docstring + justify_docstring + return_docstring


class CategoricalFormatter(object):

    def __init__(self, categorical, buf=None, length=True, na_rep='NaN',
                 footer=True):
        self.categorical = categorical
        self.buf = buf if buf is not None else StringIO(u(""))
        self.na_rep = na_rep
        self.length = length
        self.footer = footer

    def _get_footer(self):
        footer = ''

        if self.length:
            if footer:
                footer += ', '
            footer += "Length: {length}".format(length=len(self.categorical))

        level_info = self.categorical._repr_categories_info()

        # Levels are added in a newline
        if footer:
            footer += '\n'
        footer += level_info

        return compat.text_type(footer)

    def _get_formatted_values(self):
        return format_array(self.categorical.get_values(), None,
                            float_format=None, na_rep=self.na_rep)

    def to_string(self):
        categorical = self.categorical

        if len(categorical) == 0:
            if self.footer:
                return self._get_footer()
            else:
                return u('')

        fmt_values = self._get_formatted_values()

        result = [u('{i}').format(i=i) for i in fmt_values]
        result = [i.strip() for i in result]
        result = u(', ').join(result)
        result = [u('[') + result + u(']')]
        if self.footer:
            footer = self._get_footer()
            if footer:
                result.append(footer)

        return compat.text_type(u('\n').join(result))


class SeriesFormatter(object):

    def __init__(self, series, buf=None, length=True, header=True, index=True,
                 na_rep='NaN', name=False, float_format=None, dtype=True,
                 max_rows=None):
        self.series = series
        self.buf = buf if buf is not None else StringIO()
        self.name = name
        self.na_rep = na_rep
        self.header = header
        self.length = length
        self.index = index
        self.max_rows = max_rows

        if float_format is None:
            float_format = get_option("display.float_format")
        self.float_format = float_format
        self.dtype = dtype
        self.adj = _get_adjustment()

        self._chk_truncate()

    def _chk_truncate(self):
        from pandas.core.reshape.concat import concat
        max_rows = self.max_rows
        truncate_v = max_rows and (len(self.series) > max_rows)
        series = self.series
        if truncate_v:
            if max_rows == 1:
                row_num = max_rows
                series = series.iloc[:max_rows]
            else:
                row_num = max_rows // 2
                series = concat((series.iloc[:row_num],
                                 series.iloc[-row_num:]))
            self.tr_row_num = row_num
        self.tr_series = series
        self.truncate_v = truncate_v

    def _get_footer(self):
        name = self.series.name
        footer = u('')

        if getattr(self.series.index, 'freq', None) is not None:
            footer += 'Freq: {freq}'.format(freq=self.series.index.freqstr)

        if self.name is not False and name is not None:
            if footer:
                footer += ', '

            series_name = pprint_thing(name,
                                       escape_chars=('\t', '\r', '\n'))
            footer += ((u"Name: {sname}".format(sname=series_name))
                       if name is not None else "")

        if (self.length is True or
                (self.length == 'truncate' and self.truncate_v)):
            if footer:
                footer += ', '
            footer += 'Length: {length}'.format(length=len(self.series))

        if self.dtype is not False and self.dtype is not None:
            name = getattr(self.tr_series.dtype, 'name', None)
            if name:
                if footer:
                    footer += ', '
                footer += u'dtype: {typ}'.format(typ=pprint_thing(name))

        # level infos are added to the end and in a new line, like it is done
        # for Categoricals
        if is_categorical_dtype(self.tr_series.dtype):
            level_info = self.tr_series._values._repr_categories_info()
            if footer:
                footer += "\n"
            footer += level_info

        return compat.text_type(footer)

    def _get_formatted_index(self):
        index = self.tr_series.index
        is_multi = isinstance(index, MultiIndex)

        if is_multi:
            have_header = any(name for name in index.names)
            fmt_index = index.format(names=True)
        else:
            have_header = index.name is not None
            fmt_index = index.format(name=True)
        return fmt_index, have_header

    def _get_formatted_values(self):
        values_to_format = self.tr_series._formatting_values()
        return format_array(values_to_format, None,
                            float_format=self.float_format, na_rep=self.na_rep)

    def to_string(self):
        series = self.tr_series
        footer = self._get_footer()

        if len(series) == 0:
            return 'Series([], ' + footer + ')'

        fmt_index, have_header = self._get_formatted_index()
        fmt_values = self._get_formatted_values()

        if self.truncate_v:
            n_header_rows = 0
            row_num = self.tr_row_num
            width = self.adj.len(fmt_values[row_num - 1])
            if width > 3:
                dot_str = '...'
            else:
                dot_str = '..'
            # Series uses mode=center because it has single value columns
            # DataFrame uses mode=left
            dot_str = self.adj.justify([dot_str], width, mode='center')[0]
            fmt_values.insert(row_num + n_header_rows, dot_str)
            fmt_index.insert(row_num + 1, '')

        if self.index:
            result = self.adj.adjoin(3, *[fmt_index[1:], fmt_values])
        else:
            result = self.adj.adjoin(3, fmt_values).replace('\n ',
                                                            '\n').strip()

        if self.header and have_header:
            result = fmt_index[0] + '\n' + result

        if footer:
            result += '\n' + footer

        return compat.text_type(u('').join(result))


class TextAdjustment(object):

    def __init__(self):
        self.encoding = get_option("display.encoding")

    def len(self, text):
        return compat.strlen(text, encoding=self.encoding)

    def justify(self, texts, max_len, mode='right'):
        return justify(texts, max_len, mode=mode)

    def adjoin(self, space, *lists, **kwargs):
        return adjoin(space, *lists, strlen=self.len,
                      justfunc=self.justify, **kwargs)


class EastAsianTextAdjustment(TextAdjustment):

    def __init__(self):
        super(EastAsianTextAdjustment, self).__init__()
        if get_option("display.unicode.ambiguous_as_wide"):
            self.ambiguous_width = 2
        else:
            self.ambiguous_width = 1

    def len(self, text):
        return compat.east_asian_len(text, encoding=self.encoding,
                                     ambiguous_width=self.ambiguous_width)

    def justify(self, texts, max_len, mode='right'):
        # re-calculate padding space per str considering East Asian Width
        def _get_pad(t):
            return max_len - self.len(t) + len(t)

        if mode == 'left':
            return [x.ljust(_get_pad(x)) for x in texts]
        elif mode == 'center':
            return [x.center(_get_pad(x)) for x in texts]
        else:
            return [x.rjust(_get_pad(x)) for x in texts]


def _get_adjustment():
    use_east_asian_width = get_option("display.unicode.east_asian_width")
    if use_east_asian_width:
        return EastAsianTextAdjustment()
    else:
        return TextAdjustment()


class TableFormatter(object):

    is_truncated = False
    show_dimensions = None

    @property
    def should_show_dimensions(self):
        return (self.show_dimensions is True or
                (self.show_dimensions == 'truncate' and self.is_truncated))

    def _get_formatter(self, i):
        if isinstance(self.formatters, (list, tuple)):
            if is_integer(i):
                return self.formatters[i]
            else:
                return None
        else:
            if is_integer(i) and i not in self.columns:
                i = self.columns[i]
            return self.formatters.get(i, None)


class DataFrameFormatter(TableFormatter):
    """
    Render a DataFrame

    self.to_string() : console-friendly tabular output
    self.to_html()   : html table
    self.to_latex()   : LaTeX tabular environment table

    """

    __doc__ = __doc__ if __doc__ else ''
    __doc__ += common_docstring + justify_docstring + return_docstring

    def __init__(self, frame, buf=None, columns=None, col_space=None,
                 header=True, index=True, na_rep='NaN', formatters=None,
                 justify=None, float_format=None, sparsify=None,
                 index_names=True, line_width=None, max_rows=None,
                 max_cols=None, show_dimensions=False, decimal='.',
                 table_id=None, **kwds):
        self.frame = frame
        if buf is not None:
            self.buf = _expand_user(_stringify_path(buf))
        else:
            self.buf = StringIO()
        self.show_index_names = index_names

        if sparsify is None:
            sparsify = get_option("display.multi_sparse")

        self.sparsify = sparsify

        self.float_format = float_format
        self.formatters = formatters if formatters is not None else {}
        self.na_rep = na_rep
        self.decimal = decimal
        self.col_space = col_space
        self.header = header
        self.index = index
        self.line_width = line_width
        self.max_rows = max_rows
        self.max_cols = max_cols
        self.max_rows_displayed = min(max_rows or len(self.frame),
                                      len(self.frame))
        self.show_dimensions = show_dimensions
        self.table_id = table_id

        if justify is None:
            self.justify = get_option("display.colheader_justify")
        else:
            self.justify = justify

        self.kwds = kwds

        if columns is not None:
            self.columns = _ensure_index(columns)
            self.frame = self.frame[self.columns]
        else:
            self.columns = frame.columns

        self._chk_truncate()
        self.adj = _get_adjustment()

    def _chk_truncate(self):
        """
        Checks whether the frame should be truncated. If so, slices
        the frame up.
        """
        from pandas.core.reshape.concat import concat

        # Column of which first element is used to determine width of a dot col
        self.tr_size_col = -1

        # Cut the data to the information actually printed
        max_cols = self.max_cols
        max_rows = self.max_rows

        if max_cols == 0 or max_rows == 0:  # assume we are in the terminal
                                            # (why else = 0)
            (w, h) = get_terminal_size()
            self.w = w
            self.h = h
            if self.max_rows == 0:
                dot_row = 1
                prompt_row = 1
                if self.show_dimensions:
                    show_dimension_rows = 3
                n_add_rows = (self.header + dot_row + show_dimension_rows +
                              prompt_row)
                # rows available to fill with actual data
                max_rows_adj = self.h - n_add_rows
                self.max_rows_adj = max_rows_adj

            # Format only rows and columns that could potentially fit the
            # screen
            if max_cols == 0 and len(self.frame.columns) > w:
                max_cols = w
            if max_rows == 0 and len(self.frame) > h:
                max_rows = h

        if not hasattr(self, 'max_rows_adj'):
            self.max_rows_adj = max_rows
        if not hasattr(self, 'max_cols_adj'):
            self.max_cols_adj = max_cols

        max_cols_adj = self.max_cols_adj
        max_rows_adj = self.max_rows_adj

        truncate_h = max_cols_adj and (len(self.columns) > max_cols_adj)
        truncate_v = max_rows_adj and (len(self.frame) > max_rows_adj)

        frame = self.frame
        if truncate_h:
            if max_cols_adj == 0:
                col_num = len(frame.columns)
            elif max_cols_adj == 1:
                frame = frame.iloc[:, :max_cols]
                col_num = max_cols
            else:
                col_num = (max_cols_adj // 2)
                frame = concat((frame.iloc[:, :col_num],
                                frame.iloc[:, -col_num:]), axis=1)
            self.tr_col_num = col_num
        if truncate_v:
            if max_rows_adj == 0:
                row_num = len(frame)
            if max_rows_adj == 1:
                row_num = max_rows
                frame = frame.iloc[:max_rows, :]
            else:
                row_num = max_rows_adj // 2
                frame = concat((frame.iloc[:row_num, :],
                                frame.iloc[-row_num:, :]))
            self.tr_row_num = row_num

        self.tr_frame = frame
        self.truncate_h = truncate_h
        self.truncate_v = truncate_v
        self.is_truncated = self.truncate_h or self.truncate_v

    def _to_str_columns(self):
        """
        Render a DataFrame to a list of columns (as lists of strings).
        """
        frame = self.tr_frame

        # may include levels names also

        str_index = self._get_formatted_index(frame)

        if not is_list_like(self.header) and not self.header:
            stringified = []
            for i, c in enumerate(frame):
                fmt_values = self._format_col(i)
                fmt_values = _make_fixed_width(fmt_values, self.justify,
                                               minimum=(self.col_space or 0),
                                               adj=self.adj)
                stringified.append(fmt_values)
        else:
            if is_list_like(self.header):
                if len(self.header) != len(self.columns):
                    raise ValueError(('Writing {ncols} cols but got {nalias} '
                                      'aliases'
                                      .format(ncols=len(self.columns),
                                              nalias=len(self.header))))
                str_columns = [[label] for label in self.header]
            else:
                str_columns = self._get_formatted_column_labels(frame)

            stringified = []
            for i, c in enumerate(frame):
                cheader = str_columns[i]
                header_colwidth = max(self.col_space or 0,
                                      *(self.adj.len(x) for x in cheader))
                fmt_values = self._format_col(i)
                fmt_values = _make_fixed_width(fmt_values, self.justify,
                                               minimum=header_colwidth,
                                               adj=self.adj)

                max_len = max(max(self.adj.len(x) for x in fmt_values),
                              header_colwidth)
                cheader = self.adj.justify(cheader, max_len, mode=self.justify)
                stringified.append(cheader + fmt_values)

        strcols = stringified
        if self.index:
            strcols.insert(0, str_index)

        # Add ... to signal truncated
        truncate_h = self.truncate_h
        truncate_v = self.truncate_v

        if truncate_h:
            col_num = self.tr_col_num
            # infer from column header
            col_width = self.adj.len(strcols[self.tr_size_col][0])
            strcols.insert(self.tr_col_num + 1, ['...'.center(col_width)] *
                           (len(str_index)))
        if truncate_v:
            n_header_rows = len(str_index) - len(frame)
            row_num = self.tr_row_num
            for ix, col in enumerate(strcols):
                # infer from above row
                cwidth = self.adj.len(strcols[ix][row_num])
                is_dot_col = False
                if truncate_h:
                    is_dot_col = ix == col_num + 1
                if cwidth > 3 or is_dot_col:
                    my_str = '...'
                else:
                    my_str = '..'

                if ix == 0:
                    dot_mode = 'left'
                elif is_dot_col:
                    cwidth = self.adj.len(strcols[self.tr_size_col][0])
                    dot_mode = 'center'
                else:
                    dot_mode = 'right'
                dot_str = self.adj.justify([my_str], cwidth, mode=dot_mode)[0]
                strcols[ix].insert(row_num + n_header_rows, dot_str)
        return strcols

    def to_string(self):
        """
        Render a DataFrame to a console-friendly tabular output.
        """
        from pandas import Series

        frame = self.frame

        if len(frame.columns) == 0 or len(frame.index) == 0:
            info_line = (u('Empty {name}\nColumns: {col}\nIndex: {idx}')
                         .format(name=type(self.frame).__name__,
                         col=pprint_thing(frame.columns),
                         idx=pprint_thing(frame.index)))
            text = info_line
        else:

            strcols = self._to_str_columns()
            if self.line_width is None:  # no need to wrap around just print
                # the whole frame
                text = self.adj.adjoin(1, *strcols)
            elif (not isinstance(self.max_cols, int) or
                    self.max_cols > 0):  # need to wrap around
                text = self._join_multiline(*strcols)
            else:  # max_cols == 0. Try to fit frame to terminal
                text = self.adj.adjoin(1, *strcols).split('\n')
                max_len = Series(text).str.len().max()
                headers = [ele[0] for ele in strcols]
                # Size of last col determines dot col size. See
                # `self._to_str_columns
                size_tr_col = len(headers[self.tr_size_col])
                max_len += size_tr_col  # Need to make space for largest row
                # plus truncate dot col
                dif = max_len - self.w
                # '+ 1' to avoid too wide repr (GH PR #17023)
                adj_dif = dif + 1
                col_lens = Series([Series(ele).apply(len).max()
                                   for ele in strcols])
                n_cols = len(col_lens)
                counter = 0
                while adj_dif > 0 and n_cols > 1:
                    counter += 1
                    mid = int(round(n_cols / 2.))
                    mid_ix = col_lens.index[mid]
                    col_len = col_lens[mid_ix]
                    adj_dif -= (col_len + 1)  # adjoin adds one
                    col_lens = col_lens.drop(mid_ix)
                    n_cols = len(col_lens)
                max_cols_adj = n_cols - self.index  # subtract index column
                self.max_cols_adj = max_cols_adj

                # Call again _chk_truncate to cut frame appropriately
                # and then generate string representation
                self._chk_truncate()
                strcols = self._to_str_columns()
                text = self.adj.adjoin(1, *strcols)
        if not self.index:
            text = text.replace('\n ', '\n').strip()
        self.buf.writelines(text)

        if self.should_show_dimensions:
            self.buf.write("\n\n[{nrows} rows x {ncols} columns]"
                           .format(nrows=len(frame), ncols=len(frame.columns)))

    def _join_multiline(self, *strcols):
        lwidth = self.line_width
        adjoin_width = 1
        strcols = list(strcols)
        if self.index:
            idx = strcols.pop(0)
            lwidth -= np.array([self.adj.len(x)
                                for x in idx]).max() + adjoin_width

        col_widths = [np.array([self.adj.len(x) for x in col]).max() if
                      len(col) > 0 else 0 for col in strcols]
        col_bins = _binify(col_widths, lwidth)
        nbins = len(col_bins)

        if self.truncate_v:
            nrows = self.max_rows_adj + 1
        else:
            nrows = len(self.frame)

        str_lst = []
        st = 0
        for i, ed in enumerate(col_bins):
            row = strcols[st:ed]
            if self.index:
                row.insert(0, idx)
            if nbins > 1:
                if ed <= len(strcols) and i < nbins - 1:
                    row.append([' \\'] + ['  '] * (nrows - 1))
                else:
                    row.append([' '] * nrows)
            str_lst.append(self.adj.adjoin(adjoin_width, *row))
            st = ed
        return '\n\n'.join(str_lst)

    def to_latex(self, column_format=None, longtable=False, encoding=None,
                 multicolumn=False, multicolumn_format=None, multirow=False):
        """
        Render a DataFrame to a LaTeX tabular/longtable environment output.
        """

        from pandas.io.formats.latex import LatexFormatter
        latex_renderer = LatexFormatter(self, column_format=column_format,
                                        longtable=longtable,
                                        multicolumn=multicolumn,
                                        multicolumn_format=multicolumn_format,
                                        multirow=multirow)

        if encoding is None:
            encoding = 'ascii' if compat.PY2 else 'utf-8'

        if hasattr(self.buf, 'write'):
            latex_renderer.write_result(self.buf)
        elif isinstance(self.buf, compat.string_types):
            import codecs
            with codecs.open(self.buf, 'w', encoding=encoding) as f:
                latex_renderer.write_result(f)
        else:
            raise TypeError('buf is not a file name and it has no write '
                            'method')

    def _format_col(self, i):
        frame = self.tr_frame
        formatter = self._get_formatter(i)
        values_to_format = frame.iloc[:, i]._formatting_values()
        return format_array(values_to_format, formatter,
                            float_format=self.float_format, na_rep=self.na_rep,
                            space=self.col_space, decimal=self.decimal)

    def to_html(self, classes=None, notebook=False, border=None):
        """
        Render a DataFrame to a html table.

        Parameters
        ----------
        classes : str or list-like
            classes to include in the `class` attribute of the opening
            ``<table>`` tag, in addition to the default "dataframe".
        notebook : {True, False}, optional, default False
            Whether the generated HTML is for IPython Notebook.
        border : int
            A ``border=border`` attribute is included in the opening
            ``<table>`` tag. Default ``pd.options.html.border``.

            .. versionadded:: 0.19.0
         """
        from pandas.io.formats.html import HTMLFormatter
        html_renderer = HTMLFormatter(self, classes=classes,
                                      max_rows=self.max_rows,
                                      max_cols=self.max_cols,
                                      notebook=notebook,
                                      border=border,
                                      table_id=self.table_id)
        if hasattr(self.buf, 'write'):
            html_renderer.write_result(self.buf)
        elif isinstance(self.buf, compat.string_types):
            with open(self.buf, 'w') as f:
                html_renderer.write_result(f)
        else:
            raise TypeError('buf is not a file name and it has no write '
                            ' method')

    def _get_formatted_column_labels(self, frame):
        from pandas.core.index import _sparsify

        columns = frame.columns

        if isinstance(columns, MultiIndex):
            fmt_columns = columns.format(sparsify=False, adjoin=False)
            fmt_columns = lzip(*fmt_columns)
            dtypes = self.frame.dtypes._values

            # if we have a Float level, they don't use leading space at all
            restrict_formatting = any(l.is_floating for l in columns.levels)
            need_leadsp = dict(zip(fmt_columns, map(is_numeric_dtype, dtypes)))

            def space_format(x, y):
                if (y not in self.formatters and
                        need_leadsp[x] and not restrict_formatting):
                    return ' ' + y
                return y

            str_columns = list(zip(*[[space_format(x, y) for y in x]
                                     for x in fmt_columns]))
            if self.sparsify:
                str_columns = _sparsify(str_columns)

            str_columns = [list(x) for x in zip(*str_columns)]
        else:
            fmt_columns = columns.format()
            dtypes = self.frame.dtypes
            need_leadsp = dict(zip(fmt_columns, map(is_numeric_dtype, dtypes)))
            str_columns = [[' ' + x if not self._get_formatter(i) and
                            need_leadsp[x] else x]
                           for i, (col, x) in enumerate(zip(columns,
                                                            fmt_columns))]

        if self.show_index_names and self.has_index_names:
            for x in str_columns:
                x.append('')

        # self.str_columns = str_columns
        return str_columns

    @property
    def has_index_names(self):
        return _has_names(self.frame.index)

    @property
    def has_column_names(self):
        return _has_names(self.frame.columns)

    def _get_formatted_index(self, frame):
        # Note: this is only used by to_string() and to_latex(), not by
        # to_html().
        index = frame.index
        columns = frame.columns

        show_index_names = self.show_index_names and self.has_index_names
        show_col_names = (self.show_index_names and self.has_column_names)

        fmt = self._get_formatter('__index__')

        if isinstance(index, MultiIndex):
            fmt_index = index.format(sparsify=self.sparsify, adjoin=False,
                                     names=show_index_names, formatter=fmt)
        else:
            fmt_index = [index.format(name=show_index_names, formatter=fmt)]
        fmt_index = [tuple(_make_fixed_width(list(x), justify='left',
                                             minimum=(self.col_space or 0),
                                             adj=self.adj)) for x in fmt_index]

        adjoined = self.adj.adjoin(1, *fmt_index).split('\n')

        # empty space for columns
        if show_col_names:
            col_header = ['{x}'.format(x=x)
                          for x in self._get_column_name_list()]
        else:
            col_header = [''] * columns.nlevels

        if self.header:
            return col_header + adjoined
        else:
            return adjoined

    def _get_column_name_list(self):
        names = []
        columns = self.frame.columns
        if isinstance(columns, MultiIndex):
            names.extend('' if name is None else name
                         for name in columns.names)
        else:
            names.append('' if columns.name is None else columns.name)
        return names

# ----------------------------------------------------------------------
# Array formatters


def format_array(values, formatter, float_format=None, na_rep='NaN',
                 digits=None, space=None, justify='right', decimal='.'):

    if is_categorical_dtype(values):
        fmt_klass = CategoricalArrayFormatter
    elif is_interval_dtype(values):
        fmt_klass = IntervalArrayFormatter
    elif is_float_dtype(values.dtype):
        fmt_klass = FloatArrayFormatter
    elif is_period_arraylike(values):
        fmt_klass = PeriodArrayFormatter
    elif is_integer_dtype(values.dtype):
        fmt_klass = IntArrayFormatter
    elif is_datetimetz(values):
        fmt_klass = Datetime64TZFormatter
    elif is_datetime64_dtype(values.dtype):
        fmt_klass = Datetime64Formatter
    elif is_timedelta64_dtype(values.dtype):
        fmt_klass = Timedelta64Formatter
    else:
        fmt_klass = GenericArrayFormatter

    if space is None:
        space = get_option("display.column_space")

    if float_format is None:
        float_format = get_option("display.float_format")

    if digits is None:
        digits = get_option("display.precision")

    fmt_obj = fmt_klass(values, digits=digits, na_rep=na_rep,
                        float_format=float_format, formatter=formatter,
                        space=space, justify=justify, decimal=decimal)

    return fmt_obj.get_result()


class GenericArrayFormatter(object):

    def __init__(self, values, digits=7, formatter=None, na_rep='NaN',
                 space=12, float_format=None, justify='right', decimal='.',
                 quoting=None, fixed_width=True):
        self.values = values
        self.digits = digits
        self.na_rep = na_rep
        self.space = space
        self.formatter = formatter
        self.float_format = float_format
        self.justify = justify
        self.decimal = decimal
        self.quoting = quoting
        self.fixed_width = fixed_width

    def get_result(self):
        fmt_values = self._format_strings()
        return _make_fixed_width(fmt_values, self.justify)

    def _format_strings(self):
        if self.float_format is None:
            float_format = get_option("display.float_format")
            if float_format is None:
                fmt_str = ('{{x: .{prec:d}g}}'
                           .format(prec=get_option("display.precision")))
                float_format = lambda x: fmt_str.format(x=x)
        else:
            float_format = self.float_format

        formatter = (
            self.formatter if self.formatter is not None else
            (lambda x: pprint_thing(x, escape_chars=('\t', '\r', '\n'))))

        def _format(x):
            if self.na_rep is not None and is_scalar(x) and isna(x):
                if x is None:
                    return 'None'
                elif x is pd.NaT:
                    return 'NaT'
                return self.na_rep
            elif isinstance(x, PandasObject):
                return u'{x}'.format(x=x)
            else:
                # object dtype
                return u'{x}'.format(x=formatter(x))

        vals = self.values
        if isinstance(vals, Index):
            vals = vals._values
        elif isinstance(vals, ABCSparseArray):
            vals = vals.values

        is_float_type = lib.map_infer(vals, is_float) & notna(vals)
        leading_space = is_float_type.any()

        fmt_values = []
        for i, v in enumerate(vals):
            if not is_float_type[i] and leading_space:
                fmt_values.append(u' {v}'.format(v=_format(v)))
            elif is_float_type[i]:
                fmt_values.append(float_format(v))
            else:
                fmt_values.append(u' {v}'.format(v=_format(v)))

        return fmt_values


class FloatArrayFormatter(GenericArrayFormatter):
    """

    """

    def __init__(self, *args, **kwargs):
        GenericArrayFormatter.__init__(self, *args, **kwargs)

        # float_format is expected to be a string
        # formatter should be used to pass a function
        if self.float_format is not None and self.formatter is None:
            if callable(self.float_format):
                self.formatter = self.float_format
                self.float_format = None

    def _value_formatter(self, float_format=None, threshold=None):
        """Returns a function to be applied on each value to format it
        """

        # the float_format parameter supersedes self.float_format
        if float_format is None:
            float_format = self.float_format

        # we are going to compose different functions, to first convert to
        # a string, then replace the decimal symbol, and finally chop according
        # to the threshold

        # when there is no float_format, we use str instead of '%g'
        # because str(0.0) = '0.0' while '%g' % 0.0 = '0'
        if float_format:
            def base_formatter(v):
                return float_format(value=v) if notna(v) else self.na_rep
        else:
            def base_formatter(v):
                return str(v) if notna(v) else self.na_rep

        if self.decimal != '.':
            def decimal_formatter(v):
                return base_formatter(v).replace('.', self.decimal, 1)
        else:
            decimal_formatter = base_formatter

        if threshold is None:
            return decimal_formatter

        def formatter(value):
            if notna(value):
                if abs(value) > threshold:
                    return decimal_formatter(value)
                else:
                    return decimal_formatter(0.0)
            else:
                return self.na_rep

        return formatter

    def get_result_as_array(self):
        """
        Returns the float values converted into strings using
        the parameters given at initialisation, as a numpy array
        """

        if self.formatter is not None:
            return np.array([self.formatter(x) for x in self.values])

        if self.fixed_width:
            threshold = get_option("display.chop_threshold")
        else:
            threshold = None

        # if we have a fixed_width, we'll need to try different float_format
        def format_values_with(float_format):
            formatter = self._value_formatter(float_format, threshold)

            # separate the wheat from the chaff
            values = self.values
            mask = isna(values)
            if hasattr(values, 'to_dense'):  # sparse numpy ndarray
                values = values.to_dense()
            values = np.array(values, dtype='object')
            values[mask] = self.na_rep
            imask = (~mask).ravel()
            values.flat[imask] = np.array([formatter(val)
                                           for val in values.ravel()[imask]])

            if self.fixed_width:
                return _trim_zeros(values, self.na_rep)

            return values

        # There is a special default string when we are fixed-width
        # The default is otherwise to use str instead of a formatting string
        if self.float_format is None:
            if self.fixed_width:
                float_format = partial('{value: .{digits:d}f}'.format,
                                       digits=self.digits)
            else:
                float_format = self.float_format
        else:
            float_format = lambda value: self.float_format % value

        formatted_values = format_values_with(float_format)

        if not self.fixed_width:
            return formatted_values

        # we need do convert to engineering format if some values are too small
        # and would appear as 0, or if some values are too big and take too
        # much space

        if len(formatted_values) > 0:
            maxlen = max(len(x) for x in formatted_values)
            too_long = maxlen > self.digits + 6
        else:
            too_long = False

        with np.errstate(invalid='ignore'):
            abs_vals = np.abs(self.values)
            # this is pretty arbitrary for now
            # large values: more that 8 characters including decimal symbol
            # and first digit, hence > 1e6
            has_large_values = (abs_vals > 1e6).any()
            has_small_values = ((abs_vals < 10**(-self.digits)) &
                                (abs_vals > 0)).any()

        if has_small_values or (too_long and has_large_values):
            float_format = partial('{value: .{digits:d}e}'.format,
                                   digits=self.digits)
            formatted_values = format_values_with(float_format)

        return formatted_values

    def _format_strings(self):
        # shortcut
        if self.formatter is not None:
            return [self.formatter(x) for x in self.values]

        return list(self.get_result_as_array())


class IntArrayFormatter(GenericArrayFormatter):

    def _format_strings(self):
        formatter = self.formatter or (lambda x: '{x: d}'.format(x=x))
        fmt_values = [formatter(x) for x in self.values]
        return fmt_values


class Datetime64Formatter(GenericArrayFormatter):

    def __init__(self, values, nat_rep='NaT', date_format=None, **kwargs):
        super(Datetime64Formatter, self).__init__(values, **kwargs)
        self.nat_rep = nat_rep
        self.date_format = date_format

    def _format_strings(self):
        """ we by definition have DO NOT have a TZ """

        values = self.values

        if not isinstance(values, DatetimeIndex):
            values = DatetimeIndex(values)

        if self.formatter is not None and callable(self.formatter):
            return [self.formatter(x) for x in values]

        fmt_values = format_array_from_datetime(
            values.asi8.ravel(),
            format=_get_format_datetime64_from_values(values,
                                                      self.date_format),
            na_rep=self.nat_rep).reshape(values.shape)
        return fmt_values.tolist()


class IntervalArrayFormatter(GenericArrayFormatter):

    def __init__(self, values, *args, **kwargs):
        GenericArrayFormatter.__init__(self, values, *args, **kwargs)

    def _format_strings(self):
        formatter = self.formatter or str
        fmt_values = np.array([formatter(x) for x in self.values])
        return fmt_values


class PeriodArrayFormatter(IntArrayFormatter):

    def _format_strings(self):
        from pandas.core.indexes.period import IncompatibleFrequency
        try:
            values = PeriodIndex(self.values).to_native_types()
        except IncompatibleFrequency:
            # periods may contains different freq
            values = Index(self.values, dtype='object').to_native_types()

        formatter = self.formatter or (lambda x: '{x}'.format(x=x))
        fmt_values = [formatter(x) for x in values]
        return fmt_values


class CategoricalArrayFormatter(GenericArrayFormatter):

    def __init__(self, values, *args, **kwargs):
        GenericArrayFormatter.__init__(self, values, *args, **kwargs)

    def _format_strings(self):
        fmt_values = format_array(self.values.get_values(), self.formatter,
                                  float_format=self.float_format,
                                  na_rep=self.na_rep, digits=self.digits,
                                  space=self.space, justify=self.justify)
        return fmt_values


def format_percentiles(percentiles):
    """
    Outputs rounded and formatted percentiles.

    Parameters
    ----------
    percentiles : list-like, containing floats from interval [0,1]

    Returns
    -------
    formatted : list of strings

    Notes
    -----
    Rounding precision is chosen so that: (1) if any two elements of
    ``percentiles`` differ, they remain different after rounding
    (2) no entry is *rounded* to 0% or 100%.
    Any non-integer is always rounded to at least 1 decimal place.

    Examples
    --------
    Keeps all entries different after rounding:

    >>> format_percentiles([0.01999, 0.02001, 0.5, 0.666666, 0.9999])
    ['1.999%', '2.001%', '50%', '66.667%', '99.99%']

    No element is rounded to 0% or 100% (unless already equal to it).
    Duplicates are allowed:

    >>> format_percentiles([0, 0.5, 0.02001, 0.5, 0.666666, 0.9999])
    ['0%', '50%', '2.0%', '50%', '66.67%', '99.99%']
    """

    percentiles = np.asarray(percentiles)

    # It checks for np.NaN as well
    with np.errstate(invalid='ignore'):
        if not is_numeric_dtype(percentiles) or not np.all(percentiles >= 0) \
                or not np.all(percentiles <= 1):
            raise ValueError("percentiles should all be in the interval [0,1]")

    percentiles = 100 * percentiles
    int_idx = (percentiles.astype(int) == percentiles)

    if np.all(int_idx):
        out = percentiles.astype(int).astype(str)
        return [i + '%' for i in out]

    unique_pcts = np.unique(percentiles)
    to_begin = unique_pcts[0] if unique_pcts[0] > 0 else None
    to_end = 100 - unique_pcts[-1] if unique_pcts[-1] < 100 else None

    # Least precision that keeps percentiles unique after rounding
    prec = -np.floor(np.log10(np.min(
        np.ediff1d(unique_pcts, to_begin=to_begin, to_end=to_end)
    ))).astype(int)
    prec = max(1, prec)
    out = np.empty_like(percentiles, dtype=object)
    out[int_idx] = percentiles[int_idx].astype(int).astype(str)
    out[~int_idx] = percentiles[~int_idx].round(prec).astype(str)
    return [i + '%' for i in out]


def _is_dates_only(values):
    # return a boolean if we are only dates (and don't have a timezone)
    values = DatetimeIndex(values)
    if values.tz is not None:
        return False

    values_int = values.asi8
    consider_values = values_int != iNaT
    one_day_nanos = (86400 * 1e9)
    even_days = np.logical_and(consider_values,
                               values_int % int(one_day_nanos) != 0).sum() == 0
    if even_days:
        return True
    return False


def _format_datetime64(x, tz=None, nat_rep='NaT'):
    if x is None or (is_scalar(x) and isna(x)):
        return nat_rep

    if tz is not None or not isinstance(x, Timestamp):
        x = Timestamp(x, tz=tz)

    return str(x)


def _format_datetime64_dateonly(x, nat_rep='NaT', date_format=None):
    if x is None or (is_scalar(x) and isna(x)):
        return nat_rep

    if not isinstance(x, Timestamp):
        x = Timestamp(x)

    if date_format:
        return x.strftime(date_format)
    else:
        return x._date_repr


def _get_format_datetime64(is_dates_only, nat_rep='NaT', date_format=None):

    if is_dates_only:
        return lambda x, tz=None: _format_datetime64_dateonly(
            x, nat_rep=nat_rep, date_format=date_format)
    else:
        return lambda x, tz=None: _format_datetime64(x, tz=tz, nat_rep=nat_rep)


def _get_format_datetime64_from_values(values, date_format):
    """ given values and a date_format, return a string format """
    is_dates_only = _is_dates_only(values)
    if is_dates_only:
        return date_format or "%Y-%m-%d"
    return date_format


class Datetime64TZFormatter(Datetime64Formatter):

    def _format_strings(self):
        """ we by definition have a TZ """

        values = self.values.astype(object)
        is_dates_only = _is_dates_only(values)
        formatter = (self.formatter or
                     _get_format_datetime64(is_dates_only,
                                            date_format=self.date_format))
        fmt_values = [formatter(x) for x in values]

        return fmt_values


class Timedelta64Formatter(GenericArrayFormatter):

    def __init__(self, values, nat_rep='NaT', box=False, **kwargs):
        super(Timedelta64Formatter, self).__init__(values, **kwargs)
        self.nat_rep = nat_rep
        self.box = box

    def _format_strings(self):
        formatter = (self.formatter or
                     _get_format_timedelta64(self.values, nat_rep=self.nat_rep,
                                             box=self.box))
        fmt_values = np.array([formatter(x) for x in self.values])
        return fmt_values


def _get_format_timedelta64(values, nat_rep='NaT', box=False):
    """
    Return a formatter function for a range of timedeltas.
    These will all have the same format argument

    If box, then show the return in quotes
    """

    values_int = values.astype(np.int64)

    consider_values = values_int != iNaT

    one_day_nanos = (86400 * 1e9)
    even_days = np.logical_and(consider_values,
                               values_int % one_day_nanos != 0).sum() == 0
    all_sub_day = np.logical_and(
        consider_values, np.abs(values_int) >= one_day_nanos).sum() == 0

    if even_days:
        format = None
    elif all_sub_day:
        format = 'sub_day'
    else:
        format = 'long'

    def _formatter(x):
        if x is None or (is_scalar(x) and isna(x)):
            return nat_rep

        if not isinstance(x, Timedelta):
            x = Timedelta(x)
        result = x._repr_base(format=format)
        if box:
            result = "'{res}'".format(res=result)
        return result

    return _formatter


def _make_fixed_width(strings, justify='right', minimum=None, adj=None):

    if len(strings) == 0 or justify == 'all':
        return strings

    if adj is None:
        adj = _get_adjustment()

    max_len = max(adj.len(x) for x in strings)

    if minimum is not None:
        max_len = max(minimum, max_len)

    conf_max = get_option("display.max_colwidth")
    if conf_max is not None and max_len > conf_max:
        max_len = conf_max

    def just(x):
        if conf_max is not None:
            if (conf_max > 3) & (adj.len(x) > max_len):
                x = x[:max_len - 3] + '...'
        return x

    strings = [just(x) for x in strings]
    result = adj.justify(strings, max_len, mode=justify)
    return result


def _trim_zeros(str_floats, na_rep='NaN'):
    """
    Trims zeros, leaving just one before the decimal points if need be.
    """
    trimmed = str_floats

    def _cond(values):
        non_na = [x for x in values if x != na_rep]
        return (len(non_na) > 0 and all(x.endswith('0') for x in non_na) and
                not (any(('e' in x) or ('E' in x) for x in non_na)))

    while _cond(trimmed):
        trimmed = [x[:-1] if x != na_rep else x for x in trimmed]

    # leave one 0 after the decimal points if need be.
    return [x + "0" if x.endswith('.') and x != na_rep else x for x in trimmed]


def _has_names(index):
    if isinstance(index, MultiIndex):
        return com._any_not_none(*index.names)
    else:
        return index.name is not None


class EngFormatter(object):
    """
    Formats float values according to engineering format.

    Based on matplotlib.ticker.EngFormatter
    """

    # The SI engineering prefixes
    ENG_PREFIXES = {
        -24: "y",
        -21: "z",
        -18: "a",
        -15: "f",
        -12: "p",
        -9: "n",
        -6: "u",
        -3: "m",
        0: "",
        3: "k",
        6: "M",
        9: "G",
        12: "T",
        15: "P",
        18: "E",
        21: "Z",
        24: "Y"
    }

    def __init__(self, accuracy=None, use_eng_prefix=False):
        self.accuracy = accuracy
        self.use_eng_prefix = use_eng_prefix

    def __call__(self, num):
        """ Formats a number in engineering notation, appending a letter
        representing the power of 1000 of the original number. Some examples:

        >>> format_eng(0)       # for self.accuracy = 0
        ' 0'

        >>> format_eng(1000000) # for self.accuracy = 1,
                                #     self.use_eng_prefix = True
        ' 1.0M'

        >>> format_eng("-1e-6") # for self.accuracy = 2
                                #     self.use_eng_prefix = False
        '-1.00E-06'

        @param num: the value to represent
        @type num: either a numeric value or a string that can be converted to
                   a numeric value (as per decimal.Decimal constructor)

        @return: engineering formatted string
        """
        import decimal
        import math
        dnum = decimal.Decimal(str(num))

        if decimal.Decimal.is_nan(dnum):
            return 'NaN'

        if decimal.Decimal.is_infinite(dnum):
            return 'inf'

        sign = 1

        if dnum < 0:  # pragma: no cover
            sign = -1
            dnum = -dnum

        if dnum != 0:
            pow10 = decimal.Decimal(int(math.floor(dnum.log10() / 3) * 3))
        else:
            pow10 = decimal.Decimal(0)

        pow10 = pow10.min(max(self.ENG_PREFIXES.keys()))
        pow10 = pow10.max(min(self.ENG_PREFIXES.keys()))
        int_pow10 = int(pow10)

        if self.use_eng_prefix:
            prefix = self.ENG_PREFIXES[int_pow10]
        else:
            if int_pow10 < 0:
                prefix = 'E-{pow10:02d}'.format(pow10=-int_pow10)
            else:
                prefix = 'E+{pow10:02d}'.format(pow10=int_pow10)

        mant = sign * dnum / (10**pow10)

        if self.accuracy is None:  # pragma: no cover
            format_str = u("{mant: g}{prefix}")
        else:
            format_str = (u("{{mant: .{acc:d}f}}{{prefix}}")
                          .format(acc=self.accuracy))

        formatted = format_str.format(mant=mant, prefix=prefix)

        return formatted  # .strip()


def set_eng_float_format(accuracy=3, use_eng_prefix=False):
    """
    Alter default behavior on how float is formatted in DataFrame.
    Format float in engineering format. By accuracy, we mean the number of
    decimal digits after the floating point.

    See also EngFormatter.
    """

    set_option("display.float_format", EngFormatter(accuracy, use_eng_prefix))
    set_option("display.column_space", max(12, accuracy + 9))


def _binify(cols, line_width):
    adjoin_width = 1
    bins = []
    curr_width = 0
    i_last_column = len(cols) - 1
    for i, w in enumerate(cols):
        w_adjoined = w + adjoin_width
        curr_width += w_adjoined
        if i_last_column == i:
            wrap = curr_width + 1 > line_width and i > 0
        else:
            wrap = curr_width + 2 > line_width and i > 0
        if wrap:
            bins.append(i)
            curr_width = w_adjoined

    bins.append(len(cols))
    return bins


def get_level_lengths(levels, sentinel=''):
    """For each index in each level the function returns lengths of indexes.

    Parameters
    ----------
    levels : list of lists
        List of values on for level.
    sentinel : string, optional
        Value which states that no new index starts on there.

    Returns
    ----------
    Returns list of maps. For each level returns map of indexes (key is index
    in row and value is length of index).
    """
    if len(levels) == 0:
        return []

    control = [True for x in levels[0]]

    result = []
    for level in levels:
        last_index = 0

        lengths = {}
        for i, key in enumerate(level):
            if control[i] and key == sentinel:
                pass
            else:
                control[i] = False
                lengths[last_index] = i - last_index
                last_index = i

        lengths[last_index] = len(level) - last_index

        result.append(lengths)

    return result


def buffer_put_lines(buf, lines):
    """
    Appends lines to a buffer.

    Parameters
    ----------
    buf
        The buffer to write to
    lines
        The lines to append.
    """
    if any(isinstance(x, compat.text_type) for x in lines):
        lines = [compat.text_type(x) for x in lines]
    buf.write('\n'.join(lines))
