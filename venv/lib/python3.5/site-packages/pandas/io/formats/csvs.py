# -*- coding: utf-8 -*-
"""
Module for formatting output data into CSV files.
"""

from __future__ import print_function

import csv as csvlib
import numpy as np

from pandas.core.dtypes.missing import notna
from pandas.core.index import Index, MultiIndex
from pandas import compat
from pandas.compat import (StringIO, range, zip)

from pandas.io.common import (_get_handle, UnicodeWriter, _expand_user,
                              _stringify_path)
from pandas._libs import writers as libwriters
from pandas.core.indexes.datetimes import DatetimeIndex
from pandas.core.indexes.period import PeriodIndex


class CSVFormatter(object):

    def __init__(self, obj, path_or_buf=None, sep=",", na_rep='',
                 float_format=None, cols=None, header=True, index=True,
                 index_label=None, mode='w', nanRep=None, encoding=None,
                 compression=None, quoting=None, line_terminator='\n',
                 chunksize=None, tupleize_cols=False, quotechar='"',
                 date_format=None, doublequote=True, escapechar=None,
                 decimal='.'):

        self.obj = obj

        if path_or_buf is None:
            path_or_buf = StringIO()

        self.path_or_buf = _expand_user(_stringify_path(path_or_buf))
        self.sep = sep
        self.na_rep = na_rep
        self.float_format = float_format
        self.decimal = decimal

        self.header = header
        self.index = index
        self.index_label = index_label
        self.mode = mode
        self.encoding = encoding
        self.compression = compression

        if quoting is None:
            quoting = csvlib.QUOTE_MINIMAL
        self.quoting = quoting

        if quoting == csvlib.QUOTE_NONE:
            # prevents crash in _csv
            quotechar = None
        self.quotechar = quotechar

        self.doublequote = doublequote
        self.escapechar = escapechar

        self.line_terminator = line_terminator

        self.date_format = date_format

        self.tupleize_cols = tupleize_cols
        self.has_mi_columns = (isinstance(obj.columns, MultiIndex) and
                               not self.tupleize_cols)

        # validate mi options
        if self.has_mi_columns:
            if cols is not None:
                raise TypeError("cannot specify cols with a MultiIndex on the "
                                "columns")

        if cols is not None:
            if isinstance(cols, Index):
                cols = cols.to_native_types(na_rep=na_rep,
                                            float_format=float_format,
                                            date_format=date_format,
                                            quoting=self.quoting)
            else:
                cols = list(cols)
            self.obj = self.obj.loc[:, cols]

        # update columns to include possible multiplicity of dupes
        # and make sure sure cols is just a list of labels
        cols = self.obj.columns
        if isinstance(cols, Index):
            cols = cols.to_native_types(na_rep=na_rep,
                                        float_format=float_format,
                                        date_format=date_format,
                                        quoting=self.quoting)
        else:
            cols = list(cols)

        # save it
        self.cols = cols

        # preallocate data 2d list
        self.blocks = self.obj._data.blocks
        ncols = sum(b.shape[0] for b in self.blocks)
        self.data = [None] * ncols

        if chunksize is None:
            chunksize = (100000 // (len(self.cols) or 1)) or 1
        self.chunksize = int(chunksize)

        self.data_index = obj.index
        if (isinstance(self.data_index, (DatetimeIndex, PeriodIndex)) and
                date_format is not None):
            self.data_index = Index([x.strftime(date_format) if notna(x) else
                                     '' for x in self.data_index])

        self.nlevels = getattr(self.data_index, 'nlevels', 1)
        if not index:
            self.nlevels = 0

    def save(self):
        # create the writer & save
        if self.encoding is None:
            if compat.PY2:
                encoding = 'ascii'
            else:
                encoding = 'utf-8'
        else:
            encoding = self.encoding

        if hasattr(self.path_or_buf, 'write'):
            f = self.path_or_buf
            close = False
        else:
            f, handles = _get_handle(self.path_or_buf, self.mode,
                                     encoding=encoding,
                                     compression=None)
            close = True if self.compression is None else False

        try:
            writer_kwargs = dict(lineterminator=self.line_terminator,
                                 delimiter=self.sep, quoting=self.quoting,
                                 doublequote=self.doublequote,
                                 escapechar=self.escapechar,
                                 quotechar=self.quotechar)
            if encoding == 'ascii':
                self.writer = csvlib.writer(f, **writer_kwargs)
            else:
                writer_kwargs['encoding'] = encoding
                self.writer = UnicodeWriter(f, **writer_kwargs)

            self._save()

        finally:
            # GH 17778 handles compression for byte strings.
            if not close and self.compression:
                f.close()
                with open(self.path_or_buf, 'r') as f:
                    data = f.read()
                f, handles = _get_handle(self.path_or_buf, self.mode,
                                         encoding=encoding,
                                         compression=self.compression)
                f.write(data)
                close = True
            if close:
                f.close()

    def _save_header(self):

        writer = self.writer
        obj = self.obj
        index_label = self.index_label
        cols = self.cols
        has_mi_columns = self.has_mi_columns
        header = self.header
        encoded_labels = []

        has_aliases = isinstance(header, (tuple, list, np.ndarray, Index))
        if not (has_aliases or self.header):
            return
        if has_aliases:
            if len(header) != len(cols):
                raise ValueError(('Writing {ncols} cols but got {nalias} '
                                 'aliases'.format(ncols=len(cols),
                                                  nalias=len(header))))
            else:
                write_cols = header
        else:
            write_cols = cols

        if self.index:
            # should write something for index label
            if index_label is not False:
                if index_label is None:
                    if isinstance(obj.index, MultiIndex):
                        index_label = []
                        for i, name in enumerate(obj.index.names):
                            if name is None:
                                name = ''
                            index_label.append(name)
                    else:
                        index_label = obj.index.name
                        if index_label is None:
                            index_label = ['']
                        else:
                            index_label = [index_label]
                elif not isinstance(index_label,
                                    (list, tuple, np.ndarray, Index)):
                    # given a string for a DF with Index
                    index_label = [index_label]

                encoded_labels = list(index_label)
            else:
                encoded_labels = []

        if not has_mi_columns or has_aliases:
            encoded_labels += list(write_cols)
            writer.writerow(encoded_labels)
        else:
            # write out the mi
            columns = obj.columns

            # write out the names for each level, then ALL of the values for
            # each level
            for i in range(columns.nlevels):

                # we need at least 1 index column to write our col names
                col_line = []
                if self.index:

                    # name is the first column
                    col_line.append(columns.names[i])

                    if isinstance(index_label, list) and len(index_label) > 1:
                        col_line.extend([''] * (len(index_label) - 1))

                col_line.extend(columns._get_level_values(i))

                writer.writerow(col_line)

            # Write out the index line if it's not empty.
            # Otherwise, we will print out an extraneous
            # blank line between the mi and the data rows.
            if encoded_labels and set(encoded_labels) != set(['']):
                encoded_labels.extend([''] * len(columns))
                writer.writerow(encoded_labels)

    def _save(self):

        self._save_header()

        nrows = len(self.data_index)

        # write in chunksize bites
        chunksize = self.chunksize
        chunks = int(nrows / chunksize) + 1

        for i in range(chunks):
            start_i = i * chunksize
            end_i = min((i + 1) * chunksize, nrows)
            if start_i >= end_i:
                break

            self._save_chunk(start_i, end_i)

    def _save_chunk(self, start_i, end_i):

        data_index = self.data_index

        # create the data for a chunk
        slicer = slice(start_i, end_i)
        for i in range(len(self.blocks)):
            b = self.blocks[i]
            d = b.to_native_types(slicer=slicer, na_rep=self.na_rep,
                                  float_format=self.float_format,
                                  decimal=self.decimal,
                                  date_format=self.date_format,
                                  quoting=self.quoting)

            for col_loc, col in zip(b.mgr_locs, d):
                # self.data is a preallocated list
                self.data[col_loc] = col

        ix = data_index.to_native_types(slicer=slicer, na_rep=self.na_rep,
                                        float_format=self.float_format,
                                        decimal=self.decimal,
                                        date_format=self.date_format,
                                        quoting=self.quoting)

        libwriters.write_csv_rows(self.data, ix, self.nlevels,
                                  self.cols, self.writer)
