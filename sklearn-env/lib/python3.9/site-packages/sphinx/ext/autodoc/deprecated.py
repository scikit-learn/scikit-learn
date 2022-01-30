"""
    sphinx.ext.autodoc.deprecated
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    The deprecated Documenters for autodoc.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import warnings
from typing import Any

from sphinx.deprecation import RemovedInSphinx50Warning
from sphinx.ext.autodoc import (AttributeDocumenter, DataDocumenter, FunctionDocumenter,
                                MethodDocumenter)


class SingledispatchFunctionDocumenter(FunctionDocumenter):
    """
    Used to be a specialized Documenter subclass for singledispatch'ed functions.

    Retained for backwards compatibility, now does the same as the FunctionDocumenter
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        warnings.warn("%s is deprecated." % self.__class__.__name__,
                      RemovedInSphinx50Warning, stacklevel=2)
        super().__init__(*args, **kwargs)


class DataDeclarationDocumenter(DataDocumenter):
    """
    Specialized Documenter subclass for data that cannot be imported
    because they are declared without initial value (refs: PEP-526).
    """
    objtype = 'datadecl'
    directivetype = 'data'
    member_order = 60

    # must be higher than AttributeDocumenter
    priority = 11

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        warnings.warn("%s is deprecated." % self.__class__.__name__,
                      RemovedInSphinx50Warning, stacklevel=2)
        super().__init__(*args, **kwargs)


class TypeVarDocumenter(DataDocumenter):
    """
    Specialized Documenter subclass for TypeVars.
    """

    objtype = 'typevar'
    directivetype = 'data'
    priority = DataDocumenter.priority + 1  # type: ignore

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        warnings.warn("%s is deprecated." % self.__class__.__name__,
                      RemovedInSphinx50Warning, stacklevel=2)
        super().__init__(*args, **kwargs)


class SingledispatchMethodDocumenter(MethodDocumenter):
    """
    Used to be a specialized Documenter subclass for singledispatch'ed methods.

    Retained for backwards compatibility, now does the same as the MethodDocumenter
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        warnings.warn("%s is deprecated." % self.__class__.__name__,
                      RemovedInSphinx50Warning, stacklevel=2)
        super().__init__(*args, **kwargs)


class InstanceAttributeDocumenter(AttributeDocumenter):
    """
    Specialized Documenter subclass for attributes that cannot be imported
    because they are instance attributes (e.g. assigned in __init__).
    """
    objtype = 'instanceattribute'
    directivetype = 'attribute'
    member_order = 60

    # must be higher than AttributeDocumenter
    priority = 11

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        warnings.warn("%s is deprecated." % self.__class__.__name__,
                      RemovedInSphinx50Warning, stacklevel=2)
        super().__init__(*args, **kwargs)


class SlotsAttributeDocumenter(AttributeDocumenter):
    """
    Specialized Documenter subclass for attributes that cannot be imported
    because they are attributes in __slots__.
    """
    objtype = 'slotsattribute'
    directivetype = 'attribute'
    member_order = 60

    # must be higher than AttributeDocumenter
    priority = 11

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        warnings.warn("%s is deprecated." % self.__class__.__name__,
                      RemovedInSphinx50Warning, stacklevel=2)
        super().__init__(*args, **kwargs)


class GenericAliasDocumenter(DataDocumenter):
    """
    Specialized Documenter subclass for GenericAliases.
    """

    objtype = 'genericalias'
    directivetype = 'data'
    priority = DataDocumenter.priority + 1  # type: ignore

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        warnings.warn("%s is deprecated." % self.__class__.__name__,
                      RemovedInSphinx50Warning, stacklevel=2)
        super().__init__(*args, **kwargs)
