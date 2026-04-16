""" Module to handle errors in Pythran. """


class PythranInternalError(Exception):

    """ Exception raise on Incorrect internal behavior in Pythran. """


class PythranCompileError(Exception):

    """ Exception raise on when Pythran fails the compile to binary step. """


class PythranSyntaxError(SyntaxError):
    def __init__(self, msg, node=None):
        SyntaxError.__init__(self, msg)
        if node:
            self.filename = getattr(node, 'filename', None)
            self.lineno = getattr(node, 'lineno', None)
            self.offset = getattr(node, 'col_offset', None)

    def __str__(self):
        loc_info = self.lineno is not None and self.offset is not None

        if self.filename and loc_info:
            with open(self.filename) as f:
                for i in range(self.lineno - 1):
                    f.readline()  # and drop it
                extra = '{}\n{}'.format(f.readline().rstrip(),
                                        " " * (self.offset) + "^~~~ (o_0)")
        else:
            extra = None

        if loc_info:
            format_header = "{}:{}:{}"
            format_args = self.lineno, self.offset, self.args[0],
        else:
            format_header = "{}:"
            format_args = self.args[0],

        r = (format_header + " error: {}").format(
                self.filename or "<unknown>",
                *format_args)

        if extra is not None:
            r += "\n----\n"
            r += extra
            r += "\n----\n"

        return r


class PythranTypeError(PythranSyntaxError):
    "A new type to distinguish general syntax errors from typing issues"


