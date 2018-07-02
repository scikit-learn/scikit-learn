# -*- coding: utf-8 -*-
# =============================================================================
# Federal University of Rio Grande do Sul (UFRGS)
# Connectionist Artificial Intelligence Laboratory (LIAC)
# Renato de Pontes Pereira - rppereira@inf.ufrgs.br
# =============================================================================
# Copyright (c) 2011 Renato de Pontes Pereira, renato.ppontes at gmail dot com
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# =============================================================================

'''
The liac-arff module implements functions to read and write ARFF files in
Python. It was created in the Connectionist Artificial Intelligence Laboratory
(LIAC), which takes place at the Federal University of Rio Grande do Sul 
(UFRGS), in Brazil.

ARFF (Attribute-Relation File Format) is an file format specially created for
describe datasets which are commonly used for machine learning experiments and
softwares. This file format was created to be used in Weka, the best 
representative software for machine learning automated experiments.

An ARFF file can be divided into two sections: header and data. The Header 
describes the metadata of the dataset, including a general description of the 
dataset, its name and its attributes. The source below is an example of a 
header section in a XOR dataset::

    % 
    % XOR Dataset
    % 
    % Created by Renato Pereira
    %            rppereira@inf.ufrgs.br
    %            http://inf.ufrgs.br/~rppereira
    % 
    % 
    @RELATION XOR

    @ATTRIBUTE input1 REAL
    @ATTRIBUTE input2 REAL
    @ATTRIBUTE y REAL

The Data section of an ARFF file describes the observations of the dataset, in 
the case of XOR dataset::

    @DATA
    0.0,0.0,0.0
    0.0,1.0,1.0
    1.0,0.0,1.0
    1.0,1.0,0.0
    % 
    % 
    % 

Notice that several lines are starting with an ``%`` symbol, denoting a 
comment, thus, lines with ``%`` at the beginning will be ignored, except by the
description part at the beginning of the file. The declarations ``@RELATION``, 
``@ATTRIBUTE``, and ``@DATA`` are all case insensitive and obligatory.

For more information and details about the ARFF file description, consult
http://www.cs.waikato.ac.nz/~ml/weka/arff.html


ARFF Files in Python
~~~~~~~~~~~~~~~~~~~~

This module uses built-ins python objects to represent a deserialized ARFF 
file. A dictionary is used as the container of the data and metadata of ARFF,
and have the following keys:

- **description**: (OPTIONAL) a string with the description of the dataset.
- **relation**: (OBLIGATORY) a string with the name of the dataset.
- **attributes**: (OBLIGATORY) a list of attributes with the following 
  template::

    (attribute_name, attribute_type)

  the attribute_name is a string, and attribute_type must be an string
  or a list of strings.
- **data**: (OBLIGATORY) a list of data instances. Each data instance must be 
  a list with values, depending on the attributes.

The above keys must follow the case which were described, i.e., the keys are 
case sensitive. The attribute type ``attribute_type`` must be one of these 
strings (they are not case sensitive): ``NUMERIC``, ``INTEGER``, ``REAL`` or 
``STRING``. For nominal attributes, the ``atribute_type`` must be a list of 
strings.

In this format, the XOR dataset presented above can be represented as a python 
object as::

    xor_dataset = {
        'description': 'XOR Dataset',
        'relation': 'XOR',
        'attributes': [
            ('input1', 'REAL'),
            ('input2', 'REAL'),
            ('y', 'REAL'),
        ],
        'data': [
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 0.0]
        ]
    }


Features
~~~~~~~~

This module provides several features, including:

- Read and write ARFF files using python built-in structures, such dictionaries
  and lists;
- Supports `scipy.sparse.coo <http://docs.scipy
  .org/doc/scipy/reference/generated/scipy.sparse.coo_matrix.html#scipy.sparse.coo_matrix>`_
  and lists of dictionaries as used by SVMLight
- Supports the following attribute types: NUMERIC, REAL, INTEGER, STRING, and
  NOMINAL;
- Has an interface similar to other built-in modules such as ``json``, or 
  ``zipfile``;
- Supports read and write the descriptions of files;
- Supports missing values and names with spaces;
- Supports unicode values and names;
- Fully compatible with Python 2.7+ and Python 3.4+;
- Under `MIT License <http://opensource.org/licenses/MIT>`_

'''
__author__ = 'Renato de Pontes Pereira, Matthias Feurer'
__author_email__ = 'renato.ppontes@gmail.com, feurerm@informatik.uni-freiburg.de'
__version__ = '2.2.2'

import csv
import re
import sys

# CONSTANTS ===================================================================
_SIMPLE_TYPES = ['NUMERIC', 'REAL', 'INTEGER', 'STRING']

_TK_DESCRIPTION = '%'
_TK_COMMENT     = '%'
_TK_RELATION    = '@RELATION'
_TK_ATTRIBUTE   = '@ATTRIBUTE'
_TK_DATA        = '@DATA'
_TK_VALUE       = ''

_RE_RELATION     = re.compile(r'^([^\{\}%,\s]*|\".*\"|\'.*\')$', re.UNICODE)
_RE_ATTRIBUTE    = re.compile(r'^(\".*\"|\'.*\'|[^\{\}%,\s]*)\s+(.+)$', re.UNICODE)
_RE_TYPE_NOMINAL = re.compile(r'^\{\s*((\".*\"|\'.*\'|\S*)\s*,\s*)*(\".*\"|\'.*\'|\S*)\s*\}$', re.UNICODE)
_RE_ESCAPE = re.compile(r'\\\'|\\\"|\\\%|[\\"\'%]')
_RE_QUOTATION_MARKS = re.compile(r''''|"''', re.UNICODE)
_RE_REPLACE_FIRST_QUOTATION_MARK = re.compile(r'^(\'|\")')
_RE_REPLACE_LAST_QUOTATION_MARK = re.compile(r'(\'|\")$')

_ESCAPE_DCT = {
    ',': ',',
    ' ': ' ',
    "'": "\\'",
    '"': '\\"',
    '%': '\\%',
    '\\': '\\',
    '\\\'': '\\\'',
    '\\"': '\\"',
    '\\%': '\\%',
}

DENSE = 0   # Constant value representing a dense matrix
COO = 1     # Constant value representing a sparse matrix in coordinate format
LOD = 2     # Constant value representing a sparse matrix in list of
            # dictionaries format
_SUPPORTED_DATA_STRUCTURES = [DENSE, COO, LOD]

# =============================================================================

# COMPATIBILITY WITH PYTHON 3 =================================================
PY3 = sys.version_info[0] == 3
if PY3:
    unicode = str
    basestring = str
    xrange = range
# COMPABILITY WITH PYTHON 2 ===================================================
# =============================================================================
PY2 = sys.version_info[0] == 2
if PY2:
    from itertools import izip as zip

# EXCEPTIONS ==================================================================
class ArffException(Exception):
    message = None

    def __init__(self):
        self.line = -1

    def __str__(self):
        return self.message%self.line

class BadRelationFormat(ArffException):
    '''Error raised when the relation declaration is in an invalid format.'''
    message = 'Bad @RELATION format, at line %d.'

class BadAttributeFormat(ArffException):
    '''Error raised when some attribute declaration is in an invalid format.'''
    message = 'Bad @ATTRIBUTE format, at line %d.'

class BadDataFormat(ArffException):
    '''Error raised when some data instance is in an invalid format.'''
    message = 'Bad @DATA instance format, at line %d.'

class BadAttributeType(ArffException):
    '''Error raised when some invalid type is provided into the attribute 
    declaration.'''
    message = 'Bad @ATTRIBUTE type, at line %d.'

class BadNominalValue(ArffException):
    '''Error raised when a value in used in some data instance but is not 
    declared into it respective attribute declaration.'''

    def __init__(self, value):
        super(BadNominalValue, self).__init__()
        self.message = ('Data value %s not found in nominal declaration, ' % value) + 'at line %d.'

class BadNumericalValue(ArffException):
    '''Error raised when and invalid numerical value is used in some data 
    instance.'''
    message = 'Invalid numerical value, at line %d.'

class BadLayout(ArffException):
    '''Error raised when the layout of the ARFF file has something wrong.'''
    message = 'Invalid layout of the ARFF file, at line %d.'

class BadObject(ArffException):
    '''Error raised when the object representing the ARFF file has something 
    wrong.'''

    def __str__(self):
        return 'Invalid object.'

class BadObject(ArffException):
    '''Error raised when the object representing the ARFF file has something 
    wrong.'''
    def __init__(self, msg=''):
        self.msg = msg

    def __str__(self):
        return '%s'%self.msg
# =============================================================================

# INTERNAL ====================================================================
def encode_string(s):
    def replace(match):
        return _ESCAPE_DCT[match.group(0)]
    return u"'" + _RE_ESCAPE.sub(replace, s) + u"'"

class Conversor(object):
    '''Conversor is a helper used for converting ARFF types to Python types.'''

    def __init__(self, type_, values=None):
        '''Contructor.'''

        self.values = values

        if type_ == 'NUMERIC' or type_ == 'REAL':
            self._conversor = self._float
        elif type_ == 'STRING':
            self._conversor = self._string
        elif type_ == 'INTEGER':
            self._conversor = self._integer
        elif type_ == 'NOMINAL':
            self._conversor = self._nominal
        elif type_ == 'ENCODED_NOMINAL':
            self._conversor = self._encoded_nominal
            self._encoded_values = {value: i for (i, value) in enumerate(values)}
        else:
            raise BadAttributeType()

    def _float(self, value):
        '''Convert the value to float.'''
        try:
            return float(value)
        except ValueError as e:
            raise BadNumericalValue()

    def _integer(self, value):
        '''Convert the value to integer.'''
        try:
            return int(float(value))
        except ValueError as e:
            raise BadNumericalValue()

    def _string(self, value):
        '''Convert the value to string.'''
        return unicode(value)

    def _nominal(self, value):
        '''Verify the value of nominal attribute and convert it to string.'''
        if value not in self.values:
            raise BadNominalValue(value)

        return self._string(value)

    def _encoded_nominal(self, value):
        '''Perform label encoding (convert labels to integers) while reading
        the .arff file.'''
        if value not in self.values:
            raise BadNominalValue(value)

        return self._encoded_values[value]

    def __call__(self, value):
        '''Convert a ``value`` to a given type. 

        This function also verify if the value is an empty string or a missing
        value, either cases, it returns None.
        '''
        value = value.strip(' ')

        if value == u'?' or value == u'':
            return None

        value = re.sub(_RE_REPLACE_FIRST_QUOTATION_MARK, '', value)
        value = re.sub(_RE_REPLACE_LAST_QUOTATION_MARK, '', value)

        return self._conversor(value)

class Data(object):
    '''Internal helper class to allow for different matrix types without
    making the code a huge collection of if statements.'''
    def __init__(self):
        self.data = []

    def decode_data(self, s, conversors):
        values = self._get_values(s)

        if values[0][0].strip(" ") == '{':
            vdict = dict(map(lambda x: (int(x[0]), x[1]),
                             [i.strip("{").strip("}").strip(" ").split(' ') for
                              i in values]))
            values = [vdict[i] if i in vdict else unicode(0) for i in
                      xrange(len(conversors))]
        # dense lines are decoded one by one
        else:
            if len(values) != len(conversors):
                raise BadDataFormat()
        values = [conversors[i](values[i]) for i in xrange(len(values))]

        self.data.append(values)

    def _get_values(self, s):
        '''(INTERNAL) Split a line into a list of values'''
        if _RE_QUOTATION_MARKS.search(s):
            return _read_csv(s.strip(' '))
        else:
            return next(csv.reader([s.strip(' ')]))

    def encode_data(self, data, attributes):
        '''(INTERNAL) Encodes a line of data.

        Data instances follow the csv format, i.e, attribute values are
        delimited by commas. After converted from csv.

        :param data: a list of values.
        :param attributes: a list of attributes. Used to check if data is valid.
        :return: a string with the encoded data line.
        '''
        for inst in data:
            if len(inst) != len(attributes):
                raise BadObject('len(inst) = {} != len(attributes) {}'.format(
                     len(inst), len(attributes),
                ))

            new_data = []
            for value in inst:
                if value is None or value == u'' or value != value:
                    s = '?'
                else:
                    s = unicode(value)
                for escape_char in _ESCAPE_DCT:
                    if escape_char in s:
                        s = encode_string(s)
                        break
                new_data.append(s)

            yield u','.join(new_data)

class COOData(Data):
    def __init__(self):
        self.data = ([], [], [])
        self._current_num_data_points = 0

    def decode_data(self, s, conversors):
        values = self._get_values(s)

        if not values[0][0].strip(" ") == '{':
            raise BadLayout()
        elif s.replace(' ', '') == '{}':
            self._current_num_data_points += 1
            return

        vdict = dict(map(lambda x: (int(x[0]), x[1]),
                         [i.strip("{").strip("}").strip(" ").split(' ')
                         for i in values]))
        col = sorted(vdict)
        values = [conversors[key](unicode(vdict[key]))
                  for key in sorted(vdict)]
        self.data[0].extend(values)
        self.data[1].extend([self._current_num_data_points] * len(values))
        self.data[2].extend(col)

        self._current_num_data_points += 1

    def encode_data(self, data, attributes):
        num_attributes = len(attributes)
        new_data = []
        current_row = 0

        row = data.row
        col = data.col
        data = data.data

        # Check if the rows are sorted
        if not all(row[i] <= row[i + 1] for i in xrange(len(row) - 1)):
            raise ValueError("liac-arff can only output COO matrices with "
                             "sorted rows.")

        for v, col, row in zip(data, col, row):
            if row > current_row:
                # Add empty rows if necessary
                while current_row < row:
                    yield " ".join([u"{", u','.join(new_data), u"}"])
                    new_data = []
                    current_row += 1

            if col >= num_attributes:
                raise BadObject()

            if v is None or v == u'' or v != v:
                s = '?'
            else:
                s = unicode(v)
            for escape_char in _ESCAPE_DCT:
                if escape_char in s:
                    s = encode_string(s)
                    break
            new_data.append("%d %s" % (col, s))

        yield " ".join([u"{", u','.join(new_data), u"}"])

class LODData(Data):
    def __init__(self):
        self.data = []

    def decode_data(self, s, conversors):
        values = self._get_values(s)

        if not values[0][0].strip(" ") == '{':
            raise BadLayout()
        elif s.replace(' ', '') == '{}':
            self.data.append({})
            return

        vdict = dict(map(lambda x: (int(x[0]), x[1]),
                         [i.strip("{").strip("}").strip(" ").split(' ')
                          for i in values]))
        for key in vdict:
            vdict[key] = conversors[key](vdict[key])
        self.data.append(vdict)

    def encode_data(self, data, attributes):
        num_attributes = len(attributes)
        for row in data:
            new_data = []

            if len(row) > 0 and max(row) >= num_attributes:
                raise BadObject()

            for col in sorted(row):
                v = row[col]
                if v is None or v == u'' or v != v:
                    s = '?'
                else:
                    s = unicode(v)
                for escape_char in _ESCAPE_DCT:
                    if escape_char in s:
                        s = encode_string(s)
                        break
                new_data.append("%d %s" % (col, s))

            yield " ".join([u"{", u','.join(new_data), u"}"])

def _get_data_object_for_decoding(matrix_type):
    if matrix_type == DENSE:
        return Data()
    elif matrix_type == COO:
        return COOData()
    elif matrix_type == LOD:
        return LODData()
    else:
        raise ValueError("Matrix type %s not supported." % str(matrix_type))

def _get_data_object_for_encoding(matrix):
    # Probably a scipy.sparse
    if hasattr(matrix, 'format'):
        if matrix.format == 'coo':
            return COOData()
    elif isinstance(matrix[0], dict):
        return LODData()
    else:
        return Data()

def _read_csv(line):
    # TODO document
    # TODO add unit tests
    # * mixed single quotes and double quotes
    # * escaped characters!
    # * does it behave like the regular csv reader?
    values = []
    quoted = False
    i = 0
    token = ''
    quote_token = False
    comma_expected = False
    only_whitespace = True

    while i < len(line):
        line_i = line[i]
        if line_i == ',' and not quoted:
            if quote_token:
                values.append(u"%s%s%s" % (quote_token, token, quote_token))
            else:
                values.append(token)
            token = ''
            i += 1
            quote_token = False
            only_whitespace = True
            comma_expected = False
        elif comma_expected:
            if line_i in (' ', '\t', '\n', '\r'):
                i += 1
            else:
                print(line_i, i, line)
                raise BadLayout()
        # Escape character
        elif line_i == '\\':
            if len(line) == i+1:
                raise BadLayout()
            if len(line) > i+2 and line[i+1] == '\\':
                # Do not trim the escape character for escaping the escape character
                token += line[i: i + 2]
            else:
                token += line[i+1: i+2]
            i += 2
        # Quoting
        elif line_i in ("'", '"') and (not quoted or line_i == quoted):
            if only_whitespace is False:
                raise ValueError(
                    'Only whitespace allowed before quoting character at '
                    'index %d in line: %s' % (i, line))
            if quoted is False:
                token = ''
                quoted = line_i
                quote_token = line_i
            elif quoted == line_i:
                quoted = False
                comma_expected = True
            else:
                raise ValueError(
                    'Inconsistent use of single quotes and double quotes for '
                    'line at character %d: %s' % (i, line)
                )
            i += 1
        elif quoted:
            token += line_i
            i += 1
        else:
            if line_i not in (' ', '\t', '\n', '\r'):
                only_whitespace = False
            token += line_i
            i += 1
    if quoted:
        raise ValueError('Quote not closed for line: %s' % line)
    if quote_token:
        values.append(u"%s%s%s" % (quote_token, token, quote_token))
    else:
        values.append(token)
    return values

# =============================================================================

# ADVANCED INTERFACE ==========================================================
class ArffDecoder(object):
    '''An ARFF decoder.'''

    def __init__(self):
        '''Constructor.'''
        self._conversors = []
        self._current_line = 0

    def _decode_comment(self, s):
        '''(INTERNAL) Decodes a comment line.

        Comments are single line strings starting, obligatorily, with the ``%``
        character, and can have any symbol, including whitespaces or special
        characters.

        This method must receive a normalized string, i.e., a string without
        padding, including the "\r\n" characters. 

        :param s: a normalized string.
        :return: a string with the decoded comment.
        '''
        res = re.sub('^\%( )?', '', s)
        return res

    def _decode_relation(self, s):
        '''(INTERNAL) Decodes a relation line.

        The relation declaration is a line with the format ``@RELATION 
        <relation-name>``, where ``relation-name`` is a string. The string must
        start with alphabetic character and must be quoted if the name includes
        spaces, otherwise this method will raise a `BadRelationFormat` exception.

        This method must receive a normalized string, i.e., a string without
        padding, including the "\r\n" characters. 

        :param s: a normalized string.
        :return: a string with the decoded relation name.
        '''
        _, v = s.split(' ', 1)
        v = v.strip()

        if not _RE_RELATION.match(v):
            raise BadRelationFormat()

        res = unicode(v.strip('"\''))
        return res

    def _decode_attribute(self, s):
        '''(INTERNAL) Decodes an attribute line.

        The attribute is the most complex declaration in an arff file. All 
        attributes must follow the template::

             @attribute <attribute-name> <datatype>

        where ``attribute-name`` is a string, quoted if the name contains any 
        whitespace, and ``datatype`` can be:

        - Numerical attributes as ``NUMERIC``, ``INTEGER`` or ``REAL``.
        - Strings as ``STRING``.
        - Dates (NOT IMPLEMENTED).
        - Nominal attributes with format:

            {<nominal-name1>, <nominal-name2>, <nominal-name3>, ...} 

        The nominal names follow the rules for the attribute names, i.e., they
        must be quoted if the name contains whitespaces.

        This method must receive a normalized string, i.e., a string without
        padding, including the "\r\n" characters. 

        :param s: a normalized string.
        :return: a tuple (ATTRIBUTE_NAME, TYPE_OR_VALUES).
        '''
        _, v = s.split(' ', 1)
        v = v.strip()

        # Verify the general structure of declaration
        m = _RE_ATTRIBUTE.match(v)
        if not m:
            raise BadAttributeFormat()

        # Extracts the raw name and type
        name, type_ = m.groups()

        # Extracts the final name
        name = unicode(name.strip('"\''))

        # Extracts the final type
        if _RE_TYPE_NOMINAL.match(type_):
            values = _read_csv(type_.strip('{} '))
            values = [
                unicode(
                    re.sub(_RE_REPLACE_LAST_QUOTATION_MARK, '',
                        re.sub(_RE_REPLACE_FIRST_QUOTATION_MARK, '', v_.strip(' '))
                    )
                )
                for v_ in values
            ]
            type_ = values

        else:
            # If not nominal, verify the type name
            type_ = unicode(type_).upper()
            if type_ not in ['NUMERIC', 'REAL', 'INTEGER', 'STRING']:
                raise BadAttributeType()

        return (name, type_)

    def _decode(self, s, encode_nominal=False, matrix_type=DENSE):
        '''Do the job the ``encode``.'''

        # Make sure this method is idempotent
        self._current_line = 0

        # If string, convert to a list of lines
        if isinstance(s, basestring):
            s = s.strip('\r\n ').replace('\r\n', '\n').split('\n')

        # Create the return object
        obj = {
            u'description': u'',
            u'relation': u'',
            u'attributes': [],
            u'data': []
        }

        # Create the data helper object
        data = _get_data_object_for_decoding(matrix_type)

        # Read all lines
        STATE = _TK_DESCRIPTION
        for row in s:
            self._current_line += 1
            # Ignore empty lines
            row = row.strip(' \r\n')
            if not row: continue

            u_row = row.upper()

            # DESCRIPTION -----------------------------------------------------
            if u_row.startswith(_TK_DESCRIPTION) and STATE == _TK_DESCRIPTION:
                obj['description'] += self._decode_comment(row) + '\n'
            # -----------------------------------------------------------------

            # RELATION --------------------------------------------------------
            elif u_row.startswith(_TK_RELATION):
                if STATE != _TK_DESCRIPTION:
                    raise BadLayout()

                STATE = _TK_RELATION
                obj['relation'] = self._decode_relation(row)
            # -----------------------------------------------------------------

            # ATTRIBUTE -------------------------------------------------------
            elif u_row.startswith(_TK_ATTRIBUTE):
                if STATE != _TK_RELATION and STATE != _TK_ATTRIBUTE:
                    raise BadLayout()

                STATE = _TK_ATTRIBUTE

                attr = self._decode_attribute(row)
                obj['attributes'].append(attr)

                if isinstance(attr[1], (list, tuple)):
                    if encode_nominal:
                        conversor = Conversor('ENCODED_NOMINAL', attr[1])
                    else:
                        conversor = Conversor('NOMINAL', attr[1])
                else:
                    conversor = Conversor(attr[1])

                self._conversors.append(conversor)
            # -----------------------------------------------------------------

            # DATA ------------------------------------------------------------
            elif u_row.startswith(_TK_DATA):
                if STATE != _TK_ATTRIBUTE:
                    raise BadLayout()

                STATE = _TK_DATA
            # -----------------------------------------------------------------

            # COMMENT ---------------------------------------------------------
            elif u_row.startswith(_TK_COMMENT):
                pass
            # -----------------------------------------------------------------

            # DATA INSTANCES --------------------------------------------------
            elif STATE == _TK_DATA:
                data.decode_data(row, self._conversors)
            # -----------------------------------------------------------------

            # UNKNOWN INFORMATION ---------------------------------------------
            else:
                raise BadLayout()
            # -----------------------------------------------------------------

        # Alter the data object
        obj['data'] = data.data
        if obj['description'].endswith('\n'):
            obj['description'] = obj['description'][:-1]

        return obj

    def decode(self, s, encode_nominal=False, return_type=DENSE):
        '''Returns the Python representation of a given ARFF file.

        When a file object is passed as an argument, this method reads lines
        iteratively, avoiding to load unnecessary information to the memory.

        :param s: a string or file object with the ARFF file.
        :param encode_nominal: boolean, if True perform a label encoding
            while reading the .arff file.
        :param return_type: determines the data structure used to store the
            dataset. Can be one of `arff.DENSE`, `arff.COO` and `arff.LOD`.
            Consult the section on `working with sparse data`_
        '''

        try:
            return self._decode(s, encode_nominal=encode_nominal,
                                matrix_type=return_type)
        except ArffException as e:
            # print e
            e.line = self._current_line
            raise e


class ArffEncoder(object):
    '''An ARFF encoder.'''

    def _encode_comment(self, s=''):
        '''(INTERNAL) Encodes a comment line.

        Comments are single line strings starting, obligatorily, with the ``%``
        character, and can have any symbol, including whitespaces or special
        characters.

        If ``s`` is None, this method will simply return an empty comment.

        :param s: (OPTIONAL) string.
        :return: a string with the encoded comment line.
        '''
        if s:
            return u'%s %s'%(_TK_COMMENT, s)
        else:
            return u'%s' % _TK_COMMENT

    def _encode_relation(self, name):
        '''(INTERNAL) Decodes a relation line.

        The relation declaration is a line with the format ``@RELATION 
        <relation-name>``, where ``relation-name`` is a string. 

        :param name: a string.
        :return: a string with the encoded relation declaration.
        '''
        for char in ' %{},':
            if char in name:
                name = '"%s"'%name
                break

        return u'%s %s'%(_TK_RELATION, name)

    def _encode_attribute(self, name, type_):
        '''(INTERNAL) Encodes an attribute line.

        The attribute follow the template::

             @attribute <attribute-name> <datatype>

        where ``attribute-name`` is a string, and ``datatype`` can be:

        - Numerical attributes as ``NUMERIC``, ``INTEGER`` or ``REAL``.
        - Strings as ``STRING``.
        - Dates (NOT IMPLEMENTED).
        - Nominal attributes with format:

            {<nominal-name1>, <nominal-name2>, <nominal-name3>, ...} 

        This method must receive a the name of the attribute and its type, if
        the attribute type is nominal, ``type`` must be a list of values.

        :param name: a string.
        :param type_: a string or a list of string.
        :return: a string with the encoded attribute declaration.
        '''
        for char in ' %{},':
            if char in name:
                name = '"%s"'%name
                break

        if isinstance(type_, (tuple, list)):
            type_tmp = []
            for i in range(len(type_)):
                type_i = type_[i]
                for escape_char in _ESCAPE_DCT:
                    if escape_char in type_[i]:
                        type_i = encode_string(type_[i])
                        break
                type_tmp.append(u'%s' % type_i)
            type_ = u'{%s}'%(u', '.join(type_tmp))

        return u'%s %s %s'%(_TK_ATTRIBUTE, name, type_)

    def encode(self, obj):
        '''Encodes a given object to an ARFF file.

        :param obj: the object containing the ARFF information.
        :return: the ARFF file as an unicode string.
        '''
        data = [row for row in self.iter_encode(obj)]

        return u'\n'.join(data)

    def iter_encode(self, obj):
        '''The iterative version of `arff.ArffEncoder.encode`.

        This encodes iteratively a given object and return, one-by-one, the 
        lines of the ARFF file.

        :param obj: the object containing the ARFF information.
        :return: (yields) the ARFF file as unicode strings.
        '''
        # DESCRIPTION
        if obj.get('description', None):
            for row in obj['description'].split('\n'):
                yield self._encode_comment(row)

        # RELATION
        if not obj.get('relation'):
            raise BadObject('Relation name not found or with invalid value.')

        yield self._encode_relation(obj['relation'])
        yield u''

        # ATTRIBUTES
        if not obj.get('attributes'):
            raise BadObject('Attributes not found.')
            
        for attr in obj['attributes']:
            # Verify for bad object format
            if not isinstance(attr, (tuple, list)) or \
               len(attr) != 2 or \
               not isinstance(attr[0], basestring):
                raise BadObject('Invalid attribute declaration "%s"'%str(attr))

            if isinstance(attr[1], basestring):
                # Verify for invalid types
                if attr[1] not in _SIMPLE_TYPES:
                    raise BadObject('Invalid attribute type "%s"'%str(attr))

            # Verify for bad object format
            elif not isinstance(attr[1], (tuple, list)):
                raise BadObject('Invalid attribute type "%s"'%str(attr))

            yield self._encode_attribute(attr[0], attr[1])
        yield u''
        attributes = obj['attributes']

        # DATA
        yield _TK_DATA
        if 'data' in obj:
            data = _get_data_object_for_encoding(obj.get('data'))
            for line in data.encode_data(obj.get('data'), attributes):
                yield line

        yield u''

# =============================================================================

# BASIC INTERFACE =============================================================
def load(fp, encode_nominal=False, return_type=DENSE):
    '''Load a file-like object containing the ARFF document and convert it into
    a Python object. 

    :param fp: a file-like object.
    :param encode_nominal: boolean, if True perform a label encoding
        while reading the .arff file.
    :param return_type: determines the data structure used to store the
        dataset. Can be one of `arff.DENSE`, `arff.COO` and `arff.LOD`.
        Consult the section on `working with sparse data`_
    :return: a dictionary.
     '''
    decoder = ArffDecoder()
    return decoder.decode(fp, encode_nominal=encode_nominal,
                          return_type=return_type)

def loads(s, encode_nominal=False, return_type=DENSE):
    '''Convert a string instance containing the ARFF document into a Python
    object.

    :param s: a string object.
    :param encode_nominal: boolean, if True perform a label encoding
        while reading the .arff file.
    :param return_type: determines the data structure used to store the
        dataset. Can be one of `arff.DENSE`, `arff.COO` and `arff.LOD`.
        Consult the section on `working with sparse data`_
    :return: a dictionary.
    '''
    decoder = ArffDecoder()
    return decoder.decode(s, encode_nominal=encode_nominal,
                          return_type=return_type)

def dump(obj, fp):
    '''Serialize an object representing the ARFF document to a given file-like 
    object.

    :param obj: a dictionary.
    :param fp: a file-like object.
    '''
    encoder = ArffEncoder()
    generator = encoder.iter_encode(obj)

    last_row = next(generator)
    for row in generator:
        fp.write(last_row + u'\n')
        last_row = row
    fp.write(last_row)

    return fp

def dumps(obj):
    '''Serialize an object representing the ARFF document, returning a string.

    :param obj: a dictionary.
    :return: a string with the ARFF document.
    '''
    encoder = ArffEncoder()
    return encoder.encode(obj)
# =============================================================================
