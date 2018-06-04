"""
Hinting the type using docstring of class/function.

It's an irreplaceable thing if you are using Dependency Injection with passive class:
http://www.martinfowler.com/articles/injection.html

Some code extracted (or based on code) from:
https://github.com/davidhalter/jedi/blob/b489019f5bd5750051122b94cc767df47751ecb7/jedi/evaluate/docstrings.py
Thanks to @davidhalter for this utils under MIT License.

Similar solutions:

    - https://www.jetbrains.com/pycharm/help/type-hinting-in-pycharm.html
    - https://www.python.org/dev/peps/pep-0484/#type-comments
    - http://www.pydev.org/manual_adv_type_hints.html
    - https://jedi.readthedocs.org/en/latest/docs/features.html#type-hinting

Discussions:

    - https://groups.google.com/d/topic/rope-dev/JlAzmZ83K1M/discussion
    - https://groups.google.com/d/topic/rope-dev/LCFNN98vckI/discussion

"""
import re

from rope.base.oi.type_hinting import utils
from rope.base.oi.type_hinting.providers import interfaces


class ParamProvider(interfaces.IParamProvider):

    def __init__(self, docstring_parser, resolver):
        """
        :type docstring_parser: rope.base.oi.type_hinting.providers.docstrings.IParamParser
        :type resolver: rope.base.oi.type_hinting.resolvers.interfaces.IResolver
        """
        self._parse_docstring = docstring_parser
        self._resolve = resolver

    def __call__(self, pyfunc, param_name):
        """
        :type pyfunc: rope.base.pyobjectsdef.PyFunction
        :type param_name: str
        :rtype: rope.base.pyobjects.PyDefinedObject | rope.base.pyobjects.PyObject or None
        """
        type_strs = self._parse_docstring(pyfunc.get_doc(), param_name)
        if type_strs:
            return self._resolve(type_strs[0], pyfunc)


class ReturnProvider(interfaces.IReturnProvider):

    def __init__(self, docstring_parser, resolver):
        """
        :type docstring_parser: rope.base.oi.type_hinting.providers.docstrings.IReturnParser
        :type resolver: rope.base.oi.type_hinting.resolvers.interfaces.IResolver
        """
        self._parse_docstring = docstring_parser
        self._resolve = resolver

    def __call__(self, pyfunc):
        """
        :type pyfunc: rope.base.pyobjectsdef.PyFunction
        :rtype: rope.base.pyobjects.PyDefinedObject | rope.base.pyobjects.PyObject or None
        """
        type_strs = self._parse_docstring(pyfunc.get_doc())
        if type_strs:
            return self._resolve(type_strs[0], pyfunc)


class AssignmentProvider(interfaces.IAssignmentProvider):

    def __init__(self, docstring_parser, resolver):
        """
        :type docstring_parser: rope.base.oi.type_hinting.providers.docstrings.IParamParser
        :type resolver: rope.base.oi.type_hinting.resolvers.interfaces.IResolver
        """
        self._parse_docstring = docstring_parser
        self._resolve = resolver

    def __call__(self, pyname):
        """
        :type pyname: rope.base.pynamesdef.AssignedName
        :rtype: rope.base.pyobjects.PyDefinedObject | rope.base.pyobjects.PyObject or None
        """
        try:
            pyclass, attr_name = utils.get_class_with_attr_name(pyname)
        except TypeError:
            return
        else:
            type_strs = self._parse_docstring(pyclass.get_doc(), attr_name)
            if type_strs:
                return self._resolve(type_strs[0], pyclass)


class IParamParser(object):

    def __call__(self, docstring, param_name):
        """
        :type docstring: str
        :type param_name: str
        """


class IReturnParser(object):

    def __call__(self, docstring):
        """
        :type docstring: str
        """


class DocstringParamParser(IParamParser):

    DOCSTRING_PARAM_PATTERNS = [
        r'\s*:type\s+%s:\s*([^\n]+)',  # Sphinx
        r'\s*:param\s+(\w+)\s+%s:[^\n]+',  # Sphinx param with type
        r'\s*@type\s+%s:\s*([^\n]+)',  # Epydoc
    ]

    def __init__(self):
        self._strip_rst_role = RSTRoleStrip()

    def __call__(self, docstring, param_name):
        """Search `docstring` for type(-s) of `param_name`.

        >>> DocstringParamParser()(':type param: int', 'param')
        ['int']
        >>> DocstringParamParser()('@type param: int', 'param')
        ['int']
        >>> DocstringParamParser()(':type param: :class:`threading.Thread`', 'param')
        ['threading.Thread']
        >>> bool(DocstringParamParser()('no document', 'param'))
        False
        >>> DocstringParamParser()(':param int param: some description', 'param')
        ['int']
        """
        if not docstring:
            return []
        patterns = [re.compile(p % re.escape(param_name))
                    for p in self.DOCSTRING_PARAM_PATTERNS]
        for pattern in patterns:
            match = pattern.search(docstring)
            if match:
                return [self._strip_rst_role(match.group(1))]

        return []


class DocstringReturnParser(IReturnParser):

    DOCSTRING_RETURN_PATTERNS = [
        re.compile(r'\s*:rtype:\s*([^\n]+)', re.M),  # Sphinx
        re.compile(r'\s*@rtype:\s*([^\n]+)', re.M),  # Epydoc
    ]

    def __init__(self):
        self._strip_rst_role = RSTRoleStrip()

    def __call__(self, docstring):
        if not docstring:
            return []
        for p in self.DOCSTRING_RETURN_PATTERNS:
            match = p.search(docstring)
            if match:
                return [self._strip_rst_role(match.group(1))]
        return []


class RSTRoleStrip(object):

    RST_ROLE_PATTERN = re.compile(r':[^`]+:`([^`]+)`')

    def __call__(self, type_str):
        """
        Strip off the part looks like a ReST role in `type_str`.

        >>> RSTRoleStrip()(':class:`ClassName`')  # strip off :class:
        'ClassName'
        >>> RSTRoleStrip()(':py:obj:`module.Object`')  # works with domain
        'module.Object'
        >>> RSTRoleStrip()('ClassName')  # do nothing when not ReST role
        'ClassName'

        See also:
        http://sphinx-doc.org/domains.html#cross-referencing-python-objects

        """
        match = self.RST_ROLE_PATTERN.match(type_str)
        if match:
            return match.group(1)
        else:
            return type_str
