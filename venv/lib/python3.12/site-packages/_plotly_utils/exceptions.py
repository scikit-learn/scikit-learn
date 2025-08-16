class PlotlyError(Exception):
    pass


class PlotlyEmptyDataError(PlotlyError):
    pass


class PlotlyGraphObjectError(PlotlyError):
    def __init__(self, message="", path=(), notes=()):
        """
        General graph object error for validation failures.

        :param (str|unicode) message: The error message.
        :param (iterable) path: A path pointing to the error.
        :param notes: Add additional notes, but keep default exception message.

        """
        self.message = message
        self.plain_message = message  # for backwards compat
        self.path = list(path)
        self.notes = notes
        super(PlotlyGraphObjectError, self).__init__(message)

    def __str__(self):
        """This is called by Python to present the error message."""
        format_dict = {
            "message": self.message,
            "path": "[" + "][".join(repr(k) for k in self.path) + "]",
            "notes": "\n".join(self.notes),
        }
        return "{message}\n\nPath To Error: {path}\n\n{notes}".format(**format_dict)


class PlotlyDictKeyError(PlotlyGraphObjectError):
    def __init__(self, obj, path, notes=()):
        """See PlotlyGraphObjectError.__init__ for param docs."""
        format_dict = {"attribute": path[-1], "object_name": obj._name}
        message = "'{attribute}' is not allowed in '{object_name}'".format(
            **format_dict
        )
        notes = [obj.help(return_help=True)] + list(notes)
        super(PlotlyDictKeyError, self).__init__(
            message=message, path=path, notes=notes
        )


class PlotlyDictValueError(PlotlyGraphObjectError):
    def __init__(self, obj, path, notes=()):
        """See PlotlyGraphObjectError.__init__ for param docs."""
        format_dict = {"attribute": path[-1], "object_name": obj._name}
        message = "'{attribute}' has invalid value inside '{object_name}'".format(
            **format_dict
        )
        notes = [obj.help(path[-1], return_help=True)] + list(notes)
        super(PlotlyDictValueError, self).__init__(
            message=message, notes=notes, path=path
        )


class PlotlyListEntryError(PlotlyGraphObjectError):
    def __init__(self, obj, path, notes=()):
        """See PlotlyGraphObjectError.__init__ for param docs."""
        format_dict = {"index": path[-1], "object_name": obj._name}
        message = "Invalid entry found in '{object_name}' at index, '{index}'".format(
            **format_dict
        )
        notes = [obj.help(return_help=True)] + list(notes)
        super(PlotlyListEntryError, self).__init__(
            message=message, path=path, notes=notes
        )


class PlotlyDataTypeError(PlotlyGraphObjectError):
    def __init__(self, obj, path, notes=()):
        """See PlotlyGraphObjectError.__init__ for param docs."""
        format_dict = {"index": path[-1], "object_name": obj._name}
        message = "Invalid entry found in '{object_name}' at index, '{index}'".format(
            **format_dict
        )
        note = "It's invalid because it doesn't contain a valid 'type' value."
        notes = [note] + list(notes)
        super(PlotlyDataTypeError, self).__init__(
            message=message, path=path, notes=notes
        )


class PlotlyKeyError(KeyError):
    """
    KeyErrors are not printed as beautifully as other errors (this is so that
    {}[''] prints    "KeyError: ''" and not "KeyError:"). So here we use
    LookupError's __str__ to make a PlotlyKeyError object which will print nicer
    error messages for KeyErrors.
    """

    def __str__(self):
        return LookupError.__str__(self)
