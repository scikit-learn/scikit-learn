# -*- coding: utf-8 -*-
"""
Module for formatting output data in HTML.
"""

from __future__ import print_function
from distutils.version import LooseVersion

from textwrap import dedent

import pandas.core.common as com
from pandas.core.index import MultiIndex
from pandas import compat
from pandas.compat import (lzip, range, map, zip, u,
                           OrderedDict, unichr)
from pandas.core.config import get_option
from pandas.io.formats.printing import pprint_thing
from pandas.io.formats.format import (get_level_lengths,
                                      buffer_put_lines)
from pandas.io.formats.format import TableFormatter


class HTMLFormatter(TableFormatter):

    indent_delta = 2

    def __init__(self, formatter, classes=None, max_rows=None, max_cols=None,
                 notebook=False, border=None, table_id=None):
        self.fmt = formatter
        self.classes = classes

        self.frame = self.fmt.frame
        self.columns = self.fmt.tr_frame.columns
        self.elements = []
        self.bold_rows = self.fmt.kwds.get('bold_rows', False)
        self.escape = self.fmt.kwds.get('escape', True)

        self.max_rows = max_rows or len(self.fmt.frame)
        self.max_cols = max_cols or len(self.fmt.columns)
        self.show_dimensions = self.fmt.show_dimensions
        self.is_truncated = (self.max_rows < len(self.fmt.frame) or
                             self.max_cols < len(self.fmt.columns))
        self.notebook = notebook
        if border is None:
            border = get_option('display.html.border')
        self.border = border
        self.table_id = table_id

    def write(self, s, indent=0):
        rs = pprint_thing(s)
        self.elements.append(' ' * indent + rs)

    def write_th(self, s, indent=0, tags=None):
        if self.fmt.col_space is not None and self.fmt.col_space > 0:
            tags = (tags or "")
            tags += ('style="min-width: {colspace};"'
                     .format(colspace=self.fmt.col_space))

        return self._write_cell(s, kind='th', indent=indent, tags=tags)

    def write_td(self, s, indent=0, tags=None):
        return self._write_cell(s, kind='td', indent=indent, tags=tags)

    def _write_cell(self, s, kind='td', indent=0, tags=None):
        if tags is not None:
            start_tag = '<{kind} {tags}>'.format(kind=kind, tags=tags)
        else:
            start_tag = '<{kind}>'.format(kind=kind)

        if self.escape:
            # escape & first to prevent double escaping of &
            esc = OrderedDict([('&', r'&amp;'), ('<', r'&lt;'),
                               ('>', r'&gt;')])
        else:
            esc = {}
        rs = pprint_thing(s, escape_chars=esc).strip()
        self.write(u'{start}{rs}</{kind}>'
                   .format(start=start_tag, rs=rs, kind=kind), indent)

    def write_tr(self, line, indent=0, indent_delta=4, header=False,
                 align=None, tags=None, nindex_levels=0):
        if tags is None:
            tags = {}

        if align is None:
            self.write('<tr>', indent)
        else:
            self.write('<tr style="text-align: {align};">'
                       .format(align=align), indent)
        indent += indent_delta

        for i, s in enumerate(line):
            val_tag = tags.get(i, None)
            if header or (self.bold_rows and i < nindex_levels):
                self.write_th(s, indent, tags=val_tag)
            else:
                self.write_td(s, indent, tags=val_tag)

        indent -= indent_delta
        self.write('</tr>', indent)

    def write_style(self):
        # We use the "scoped" attribute here so that the desired
        # style properties for the data frame are not then applied
        # throughout the entire notebook.
        template_first = """\
            <style scoped>"""
        template_last = """\
            </style>"""
        template_select = """\
                .dataframe %s {
                    %s: %s;
                }"""
        element_props = [('tbody tr th:only-of-type',
                          'vertical-align',
                          'middle'),
                         ('tbody tr th',
                          'vertical-align',
                          'top')]
        if isinstance(self.columns, MultiIndex):
            element_props.append(('thead tr th',
                                  'text-align',
                                  'left'))
            if all((self.fmt.has_index_names,
                    self.fmt.index,
                    self.fmt.show_index_names)):
                element_props.append(('thead tr:last-of-type th',
                                      'text-align',
                                      'right'))
        else:
            element_props.append(('thead th',
                                  'text-align',
                                  'right'))
        template_mid = '\n\n'.join(map(lambda t: template_select % t,
                                       element_props))
        template = dedent('\n'.join((template_first,
                                     template_mid,
                                     template_last)))
        if self.notebook:
            self.write(template)

    def write_result(self, buf):
        indent = 0
        id_section = ""
        frame = self.frame

        _classes = ['dataframe']  # Default class.
        use_mathjax = get_option("display.html.use_mathjax")
        if not use_mathjax:
            _classes.append('tex2jax_ignore')
        if self.classes is not None:
            if isinstance(self.classes, str):
                self.classes = self.classes.split()
            if not isinstance(self.classes, (list, tuple)):
                raise AssertionError('classes must be list or tuple, not {typ}'
                                     .format(typ=type(self.classes)))
            _classes.extend(self.classes)

        if self.notebook:
            div_style = ''
            try:
                import IPython
                if IPython.__version__ < LooseVersion('3.0.0'):
                    div_style = ' style="max-width:1500px;overflow:auto;"'
            except (ImportError, AttributeError):
                pass

            self.write('<div{style}>'.format(style=div_style))

        self.write_style()

        if self.table_id is not None:
            id_section = ' id="{table_id}"'.format(table_id=self.table_id)
        self.write('<table border="{border}" class="{cls}"{id_section}>'
                   .format(border=self.border, cls=' '.join(_classes),
                           id_section=id_section), indent)

        indent += self.indent_delta
        indent = self._write_header(indent)
        indent = self._write_body(indent)

        self.write('</table>', indent)
        if self.should_show_dimensions:
            by = chr(215) if compat.PY3 else unichr(215)  # ×
            self.write(u('<p>{rows} rows {by} {cols} columns</p>')
                       .format(rows=len(frame),
                               by=by,
                               cols=len(frame.columns)))

        if self.notebook:
            self.write('</div>')

        buffer_put_lines(buf, self.elements)

    def _write_header(self, indent):
        truncate_h = self.fmt.truncate_h
        row_levels = self.frame.index.nlevels
        if not self.fmt.header:
            # write nothing
            return indent

        def _column_header():
            if self.fmt.index:
                row = [''] * (self.frame.index.nlevels - 1)
            else:
                row = []

            if isinstance(self.columns, MultiIndex):
                if self.fmt.has_column_names and self.fmt.index:
                    row.append(single_column_table(self.columns.names))
                else:
                    row.append('')
                style = "text-align: {just};".format(just=self.fmt.justify)
                row.extend([single_column_table(c, self.fmt.justify, style)
                            for c in self.columns])
            else:
                if self.fmt.index:
                    row.append(self.columns.name or '')
                row.extend(self.columns)
            return row

        self.write('<thead>', indent)
        row = []

        indent += self.indent_delta

        if isinstance(self.columns, MultiIndex):
            template = 'colspan="{span:d}" halign="left"'

            if self.fmt.sparsify:
                # GH3547
                sentinel = com.sentinel_factory()
            else:
                sentinel = None
            levels = self.columns.format(sparsify=sentinel, adjoin=False,
                                         names=False)
            level_lengths = get_level_lengths(levels, sentinel)
            inner_lvl = len(level_lengths) - 1
            for lnum, (records, values) in enumerate(zip(level_lengths,
                                                         levels)):
                if truncate_h:
                    # modify the header lines
                    ins_col = self.fmt.tr_col_num
                    if self.fmt.sparsify:
                        recs_new = {}
                        # Increment tags after ... col.
                        for tag, span in list(records.items()):
                            if tag >= ins_col:
                                recs_new[tag + 1] = span
                            elif tag + span > ins_col:
                                recs_new[tag] = span + 1
                                if lnum == inner_lvl:
                                    values = (values[:ins_col] + (u('...'),) +
                                              values[ins_col:])
                                else:
                                    # sparse col headers do not receive a ...
                                    values = (values[:ins_col] +
                                              (values[ins_col - 1], ) +
                                              values[ins_col:])
                            else:
                                recs_new[tag] = span
                            # if ins_col lies between tags, all col headers
                            # get ...
                            if tag + span == ins_col:
                                recs_new[ins_col] = 1
                                values = (values[:ins_col] + (u('...'),) +
                                          values[ins_col:])
                        records = recs_new
                        inner_lvl = len(level_lengths) - 1
                        if lnum == inner_lvl:
                            records[ins_col] = 1
                    else:
                        recs_new = {}
                        for tag, span in list(records.items()):
                            if tag >= ins_col:
                                recs_new[tag + 1] = span
                            else:
                                recs_new[tag] = span
                        recs_new[ins_col] = 1
                        records = recs_new
                        values = (values[:ins_col] + [u('...')] +
                                  values[ins_col:])

                name = self.columns.names[lnum]
                row = [''] * (row_levels - 1) + ['' if name is None else
                                                 pprint_thing(name)]

                if row == [""] and self.fmt.index is False:
                    row = []

                tags = {}
                j = len(row)
                for i, v in enumerate(values):
                    if i in records:
                        if records[i] > 1:
                            tags[j] = template.format(span=records[i])
                    else:
                        continue
                    j += 1
                    row.append(v)
                self.write_tr(row, indent, self.indent_delta, tags=tags,
                              header=True)
        else:
            col_row = _column_header()
            align = self.fmt.justify

            if truncate_h:
                ins_col = row_levels + self.fmt.tr_col_num
                col_row.insert(ins_col, '...')

            self.write_tr(col_row, indent, self.indent_delta, header=True,
                          align=align)

        if all((self.fmt.has_index_names,
                self.fmt.index,
                self.fmt.show_index_names)):
            row = ([x if x is not None else ''
                    for x in self.frame.index.names] +
                   [''] * min(len(self.columns), self.max_cols))
            if truncate_h:
                ins_col = row_levels + self.fmt.tr_col_num
                row.insert(ins_col, '')
            self.write_tr(row, indent, self.indent_delta, header=True)

        indent -= self.indent_delta
        self.write('</thead>', indent)

        return indent

    def _write_body(self, indent):
        self.write('<tbody>', indent)
        indent += self.indent_delta

        fmt_values = {}
        for i in range(min(len(self.columns), self.max_cols)):
            fmt_values[i] = self.fmt._format_col(i)

        # write values
        if self.fmt.index:
            if isinstance(self.frame.index, MultiIndex):
                self._write_hierarchical_rows(fmt_values, indent)
            else:
                self._write_regular_rows(fmt_values, indent)
        else:
            for i in range(min(len(self.frame), self.max_rows)):
                row = [fmt_values[j][i] for j in range(len(self.columns))]
                self.write_tr(row, indent, self.indent_delta, tags=None)

        indent -= self.indent_delta
        self.write('</tbody>', indent)
        indent -= self.indent_delta

        return indent

    def _write_regular_rows(self, fmt_values, indent):
        truncate_h = self.fmt.truncate_h
        truncate_v = self.fmt.truncate_v

        ncols = len(self.fmt.tr_frame.columns)
        nrows = len(self.fmt.tr_frame)
        fmt = self.fmt._get_formatter('__index__')
        if fmt is not None:
            index_values = self.fmt.tr_frame.index.map(fmt)
        else:
            index_values = self.fmt.tr_frame.index.format()

        row = []
        for i in range(nrows):

            if truncate_v and i == (self.fmt.tr_row_num):
                str_sep_row = ['...' for ele in row]
                self.write_tr(str_sep_row, indent, self.indent_delta,
                              tags=None, nindex_levels=1)

            row = []
            row.append(index_values[i])
            row.extend(fmt_values[j][i] for j in range(ncols))

            if truncate_h:
                dot_col_ix = self.fmt.tr_col_num + 1
                row.insert(dot_col_ix, '...')
            self.write_tr(row, indent, self.indent_delta, tags=None,
                          nindex_levels=1)

    def _write_hierarchical_rows(self, fmt_values, indent):
        template = 'rowspan="{span}" valign="top"'

        truncate_h = self.fmt.truncate_h
        truncate_v = self.fmt.truncate_v
        frame = self.fmt.tr_frame
        ncols = len(frame.columns)
        nrows = len(frame)
        row_levels = self.frame.index.nlevels

        idx_values = frame.index.format(sparsify=False, adjoin=False,
                                        names=False)
        idx_values = lzip(*idx_values)

        if self.fmt.sparsify:
            # GH3547
            sentinel = com.sentinel_factory()
            levels = frame.index.format(sparsify=sentinel, adjoin=False,
                                        names=False)

            level_lengths = get_level_lengths(levels, sentinel)
            inner_lvl = len(level_lengths) - 1
            if truncate_v:
                # Insert ... row and adjust idx_values and
                # level_lengths to take this into account.
                ins_row = self.fmt.tr_row_num
                inserted = False
                for lnum, records in enumerate(level_lengths):
                    rec_new = {}
                    for tag, span in list(records.items()):
                        if tag >= ins_row:
                            rec_new[tag + 1] = span
                        elif tag + span > ins_row:
                            rec_new[tag] = span + 1

                            # GH 14882 - Make sure insertion done once
                            if not inserted:
                                dot_row = list(idx_values[ins_row - 1])
                                dot_row[-1] = u('...')
                                idx_values.insert(ins_row, tuple(dot_row))
                                inserted = True
                            else:
                                dot_row = list(idx_values[ins_row])
                                dot_row[inner_lvl - lnum] = u('...')
                                idx_values[ins_row] = tuple(dot_row)
                        else:
                            rec_new[tag] = span
                        # If ins_row lies between tags, all cols idx cols
                        # receive ...
                        if tag + span == ins_row:
                            rec_new[ins_row] = 1
                            if lnum == 0:
                                idx_values.insert(ins_row, tuple(
                                    [u('...')] * len(level_lengths)))

                            # GH 14882 - Place ... in correct level
                            elif inserted:
                                dot_row = list(idx_values[ins_row])
                                dot_row[inner_lvl - lnum] = u('...')
                                idx_values[ins_row] = tuple(dot_row)
                    level_lengths[lnum] = rec_new

                level_lengths[inner_lvl][ins_row] = 1
                for ix_col in range(len(fmt_values)):
                    fmt_values[ix_col].insert(ins_row, '...')
                nrows += 1

            for i in range(nrows):
                row = []
                tags = {}

                sparse_offset = 0
                j = 0
                for records, v in zip(level_lengths, idx_values[i]):
                    if i in records:
                        if records[i] > 1:
                            tags[j] = template.format(span=records[i])
                    else:
                        sparse_offset += 1
                        continue

                    j += 1
                    row.append(v)

                row.extend(fmt_values[j][i] for j in range(ncols))
                if truncate_h:
                    row.insert(row_levels - sparse_offset +
                               self.fmt.tr_col_num, '...')
                self.write_tr(row, indent, self.indent_delta, tags=tags,
                              nindex_levels=len(levels) - sparse_offset)
        else:
            for i in range(len(frame)):
                idx_values = list(zip(*frame.index.format(
                    sparsify=False, adjoin=False, names=False)))
                row = []
                row.extend(idx_values[i])
                row.extend(fmt_values[j][i] for j in range(ncols))
                if truncate_h:
                    row.insert(row_levels + self.fmt.tr_col_num, '...')
                self.write_tr(row, indent, self.indent_delta, tags=None,
                              nindex_levels=frame.index.nlevels)


def single_column_table(column, align=None, style=None):
    table = '<table'
    if align is not None:
        table += (' align="{align}"'.format(align=align))
    if style is not None:
        table += (' style="{style}"'.format(style=style))
    table += '><tbody>'
    for i in column:
        table += ('<tr><td>{i!s}</td></tr>'.format(i=i))
    table += '</tbody></table>'
    return table


def single_row_table(row):  # pragma: no cover
    table = '<table><tbody><tr>'
    for i in row:
        table += ('<td>{i!s}</td>'.format(i=i))
    table += '</tr></tbody></table>'
    return table
