import hashlib
import os

from parso._compatibility import FileNotFoundError, is_pypy
from parso.pgen2.pgen import generate_grammar
from parso.utils import split_lines, python_bytes_to_unicode, parse_version_string
from parso.python.diff import DiffParser
from parso.python.tokenize import tokenize_lines, tokenize
from parso.python import token
from parso.cache import parser_cache, load_module, save_module
from parso.parser import BaseParser
from parso.python.parser import Parser as PythonParser
from parso.python.errors import ErrorFinderConfig
from parso.python import pep8

_loaded_grammars = {}


class Grammar(object):
    """
    :py:func:`parso.load_grammar` returns instances of this class.

    Creating custom grammars by calling this is not supported, yet.
    """
    #:param text: A BNF representation of your grammar.
    _error_normalizer_config = None
    _token_namespace = None
    _default_normalizer_config = pep8.PEP8NormalizerConfig()

    def __init__(self, text, tokenizer, parser=BaseParser, diff_parser=None):
        self._pgen_grammar = generate_grammar(
            text,
            token_namespace=self._get_token_namespace()
        )
        self._parser = parser
        self._tokenizer = tokenizer
        self._diff_parser = diff_parser
        self._hashed = hashlib.sha256(text.encode("utf-8")).hexdigest()

    def parse(self, code=None, **kwargs):
        """
        If you want to parse a Python file you want to start here, most likely.

        If you need finer grained control over the parsed instance, there will be
        other ways to access it.

        :param str code: A unicode or bytes string. When it's not possible to
            decode bytes to a string, returns a
            :py:class:`UnicodeDecodeError`.
        :param bool error_recovery: If enabled, any code will be returned. If
            it is invalid, it will be returned as an error node. If disabled,
            you will get a ParseError when encountering syntax errors in your
            code.
        :param str start_symbol: The grammar symbol that you want to parse. Only
            allowed to be used when error_recovery is False.
        :param str path: The path to the file you want to open. Only needed for caching.
        :param bool cache: Keeps a copy of the parser tree in RAM and on disk
            if a path is given. Returns the cached trees if the corresponding
            files on disk have not changed.
        :param bool diff_cache: Diffs the cached python module against the new
            code and tries to parse only the parts that have changed. Returns
            the same (changed) module that is found in cache. Using this option
            requires you to not do anything anymore with the cached modules
            under that path, because the contents of it might change. This
            option is still somewhat experimental. If you want stability,
            please don't use it.
        :param bool cache_path: If given saves the parso cache in this
            directory. If not given, defaults to the default cache places on
            each platform.

        :return: A subclass of :py:class:`parso.tree.NodeOrLeaf`. Typically a
            :py:class:`parso.python.tree.Module`.
        """
        if 'start_pos' in kwargs:
            raise TypeError("parse() got an unexpected keyword argument.")
        return self._parse(code=code, **kwargs)

    def _parse(self, code=None, error_recovery=True, path=None,
               start_symbol=None, cache=False, diff_cache=False,
               cache_path=None, start_pos=(1, 0)):
        """
        Wanted python3.5 * operator and keyword only arguments. Therefore just
        wrap it all.
        start_pos here is just a parameter internally used. Might be public
        sometime in the future.
        """
        if code is None and path is None:
            raise TypeError("Please provide either code or a path.")

        if start_symbol is None:
            start_symbol = self._start_symbol

        if error_recovery and start_symbol != 'file_input':
            raise NotImplementedError("This is currently not implemented.")

        if cache and path is not None:
            module_node = load_module(self._hashed, path, cache_path=cache_path)
            if module_node is not None:
                return module_node

        if code is None:
            with open(path, 'rb') as f:
                code = f.read()

        code = python_bytes_to_unicode(code)

        lines = split_lines(code, keepends=True)
        if diff_cache:
            if self._diff_parser is None:
                raise TypeError("You have to define a diff parser to be able "
                                "to use this option.")
            try:
                module_cache_item = parser_cache[self._hashed][path]
            except KeyError:
                pass
            else:
                module_node = module_cache_item.node
                old_lines = module_cache_item.lines
                if old_lines == lines:
                    return module_node

                new_node = self._diff_parser(
                    self._pgen_grammar, self._tokenizer, module_node
                ).update(
                    old_lines=old_lines,
                    new_lines=lines
                )
                save_module(self._hashed, path, new_node, lines,
                            # Never pickle in pypy, it's slow as hell.
                            pickling=cache and not is_pypy,
                            cache_path=cache_path)
                return new_node

        tokens = self._tokenizer(lines, start_pos)

        p = self._parser(
            self._pgen_grammar,
            error_recovery=error_recovery,
            start_symbol=start_symbol
        )
        root_node = p.parse(tokens=tokens)

        if cache or diff_cache:
            save_module(self._hashed, path, root_node, lines,
                        # Never pickle in pypy, it's slow as hell.
                        pickling=cache and not is_pypy,
                        cache_path=cache_path)
        return root_node

    def _get_token_namespace(self):
        ns = self._token_namespace
        if ns is None:
            raise ValueError("The token namespace should be set.")
        return ns

    def iter_errors(self, node):
        """
        Given a :py:class:`parso.tree.NodeOrLeaf` returns a generator of
        :py:class:`parso.normalizer.Issue` objects. For Python this is
        a list of syntax/indentation errors.
        """
        if self._error_normalizer_config is None:
            raise ValueError("No error normalizer specified for this grammar.")

        return self._get_normalizer_issues(node, self._error_normalizer_config)

    def _get_normalizer(self, normalizer_config):
        if normalizer_config is None:
            normalizer_config = self._default_normalizer_config
            if normalizer_config is None:
                raise ValueError("You need to specify a normalizer, because "
                                 "there's no default normalizer for this tree.")
        return normalizer_config.create_normalizer(self)

    def _normalize(self, node, normalizer_config=None):
        """
        TODO this is not public, yet.
        The returned code will be normalized, e.g. PEP8 for Python.
        """
        normalizer = self._get_normalizer(normalizer_config)
        return normalizer.walk(node)

    def _get_normalizer_issues(self, node, normalizer_config=None):
        normalizer = self._get_normalizer(normalizer_config)
        normalizer.walk(node)
        return normalizer.issues

    def __repr__(self):
        labels = self._pgen_grammar.number2symbol.values()
        txt = ' '.join(list(labels)[:3]) + ' ...'
        return '<%s:%s>' % (self.__class__.__name__, txt)


class PythonGrammar(Grammar):
    _error_normalizer_config = ErrorFinderConfig()
    _token_namespace = token
    _start_symbol = 'file_input'

    def __init__(self, version_info, bnf_text):
        super(PythonGrammar, self).__init__(
            bnf_text,
            tokenizer=self._tokenize_lines,
            parser=PythonParser,
            diff_parser=DiffParser
        )
        self.version_info = version_info

    def _tokenize_lines(self, lines, start_pos):
        return tokenize_lines(lines, self.version_info, start_pos=start_pos)

    def _tokenize(self, code):
        # Used by Jedi.
        return tokenize(code, self.version_info)


def load_grammar(**kwargs):
    """
    Loads a :py:class:`parso.Grammar`. The default version is the current Python
    version.

    :param str version: A python version string, e.g. ``version='3.3'``.
    """
    def load_grammar(language='python', version=None):
        if language == 'python':
            version_info = parse_version_string(version)

            file = os.path.join(
                'python',
                'grammar%s%s.txt' % (version_info.major, version_info.minor)
            )

            global _loaded_grammars
            path = os.path.join(os.path.dirname(__file__), file)
            try:
                return _loaded_grammars[path]
            except KeyError:
                try:
                    with open(path) as f:
                        bnf_text = f.read()

                    grammar = PythonGrammar(version_info, bnf_text)
                    return _loaded_grammars.setdefault(path, grammar)
                except FileNotFoundError:
                    message = "Python version %s is currently not supported." % version
                    raise NotImplementedError(message)
        else:
            raise NotImplementedError("No support for language %s." % language)

    return load_grammar(**kwargs)
