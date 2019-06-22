# coding: utf-8

from __future__ import absolute_import, unicode_literals, print_function

import sys
import os
import warnings
import glob
from importlib import import_module


import ruamel.yaml
from ruamel.yaml.error import UnsafeLoaderWarning, YAMLError  # NOQA

from ruamel.yaml.tokens import *  # NOQA
from ruamel.yaml.events import *  # NOQA
from ruamel.yaml.nodes import *  # NOQA

from ruamel.yaml.loader import BaseLoader, SafeLoader, Loader, RoundTripLoader  # NOQA
from ruamel.yaml.dumper import BaseDumper, SafeDumper, Dumper, RoundTripDumper  # NOQA
from ruamel.yaml.compat import StringIO, BytesIO, with_metaclass, PY3, nprint
from ruamel.yaml.resolver import VersionedResolver, Resolver  # NOQA
from ruamel.yaml.representer import (
    BaseRepresenter,
    SafeRepresenter,
    Representer,
    RoundTripRepresenter,
)
from ruamel.yaml.constructor import (
    BaseConstructor,
    SafeConstructor,
    Constructor,
    RoundTripConstructor,
)
from ruamel.yaml.loader import Loader as UnsafeLoader

if False:  # MYPY
    from typing import List, Set, Dict, Union, Any, Callable  # NOQA
    from ruamel.yaml.compat import StreamType, StreamTextType, VersionType  # NOQA

    if PY3:
        from pathlib import Path
    else:
        Path = Any

try:
    from _ruamel_yaml import CParser, CEmitter  # type: ignore
except:  # NOQA
    CParser = CEmitter = None

# import io

enforce = object()


# YAML is an acronym, i.e. spoken: rhymes with "camel". And thus a
# subset of abbreviations, which should be all caps according to PEP8


class YAML(object):
    def __init__(
        self, _kw=enforce, typ=None, pure=False, output=None, plug_ins=None  # input=None,
    ):
        # type: (Any, Optional[Text], Any, Any, Any) -> None
        """
        _kw: not used, forces keyword arguments in 2.7 (in 3 you can do (*, safe_load=..)
        typ: 'rt'/None -> RoundTripLoader/RoundTripDumper,  (default)
             'safe'    -> SafeLoader/SafeDumper,
             'unsafe'  -> normal/unsafe Loader/Dumper
             'base'    -> baseloader
        pure: if True only use Python modules
        input/output: needed to work as context manager
        plug_ins: a list of plug-in files
        """
        if _kw is not enforce:
            raise TypeError(
                '{}.__init__() takes no positional argument but at least '
                'one was given ({!r})'.format(self.__class__.__name__, _kw)
            )

        self.typ = 'rt' if typ is None else typ
        self.pure = pure

        # self._input = input
        self._output = output
        self._context_manager = None  # type: Any

        self.plug_ins = []  # type: List[Any]
        for pu in ([] if plug_ins is None else plug_ins) + self.official_plug_ins():
            file_name = pu.replace(os.sep, '.')
            self.plug_ins.append(import_module(file_name))
        self.Resolver = ruamel.yaml.resolver.VersionedResolver  # type: Any
        self.allow_unicode = True
        self.Reader = None  # type: Any
        self.Scanner = None  # type: Any
        self.Serializer = None  # type: Any
        self.default_flow_style = None  # type: Any
        if self.typ == 'rt':
            self.default_flow_style = False
            # no optimized rt-dumper yet
            self.Emitter = ruamel.yaml.emitter.Emitter  # type: Any
            self.Serializer = ruamel.yaml.serializer.Serializer  # type: Any
            self.Representer = ruamel.yaml.representer.RoundTripRepresenter  # type: Any
            self.Scanner = ruamel.yaml.scanner.RoundTripScanner  # type: Any
            # no optimized rt-parser yet
            self.Parser = ruamel.yaml.parser.RoundTripParser  # type: Any
            self.Composer = ruamel.yaml.composer.Composer  # type: Any
            self.Constructor = ruamel.yaml.constructor.RoundTripConstructor  # type: Any
        elif self.typ == 'safe':
            self.Emitter = (
                ruamel.yaml.emitter.Emitter if pure or CEmitter is None else CEmitter
            )
            self.Representer = ruamel.yaml.representer.SafeRepresenter
            self.Parser = ruamel.yaml.parser.Parser if pure or CParser is None else CParser
            self.Composer = ruamel.yaml.composer.Composer
            self.Constructor = ruamel.yaml.constructor.SafeConstructor
        elif self.typ == 'base':
            self.Emitter = ruamel.yaml.emitter.Emitter
            self.Representer = ruamel.yaml.representer.BaseRepresenter
            self.Parser = ruamel.yaml.parser.Parser if pure or CParser is None else CParser
            self.Composer = ruamel.yaml.composer.Composer
            self.Constructor = ruamel.yaml.constructor.BaseConstructor
        elif self.typ == 'unsafe':
            self.Emitter = (
                ruamel.yaml.emitter.Emitter if pure or CEmitter is None else CEmitter
            )
            self.Representer = ruamel.yaml.representer.Representer
            self.Parser = ruamel.yaml.parser.Parser if pure or CParser is None else CParser
            self.Composer = ruamel.yaml.composer.Composer
            self.Constructor = ruamel.yaml.constructor.Constructor
        else:
            for module in self.plug_ins:
                if getattr(module, 'typ', None) == self.typ:
                    module.init_typ(self)
                    break
            else:
                raise NotImplementedError(
                    'typ "{}"not recognised (need to install plug-in?)'.format(self.typ)
                )
        self.stream = None
        self.canonical = None
        self.old_indent = None
        self.width = None
        self.line_break = None

        self.map_indent = None
        self.sequence_indent = None
        self.sequence_dash_offset = 0
        self.compact_seq_seq = None
        self.compact_seq_map = None
        self.sort_base_mapping_type_on_output = None  # default: sort

        self.top_level_colon_align = None
        self.prefix_colon = None
        self.version = None
        self.preserve_quotes = None
        self.allow_duplicate_keys = False  # duplicate keys in map, set
        self.encoding = 'utf-8'
        self.explicit_start = None
        self.explicit_end = None
        self.tags = None
        self.default_style = None
        self.top_level_block_style_scalar_no_indent_error_1_1 = False
        # [a, b: 1, c: {d: 2}]  vs. [a, {b: 1}, {c: {d: 2}}]
        self.brace_single_entry_mapping_in_flow_sequence = False

    @property
    def reader(self):
        # type: () -> Any
        try:
            return self._reader  # type: ignore
        except AttributeError:
            self._reader = self.Reader(None, loader=self)
            return self._reader

    @property
    def scanner(self):
        # type: () -> Any
        try:
            return self._scanner  # type: ignore
        except AttributeError:
            self._scanner = self.Scanner(loader=self)
            return self._scanner

    @property
    def parser(self):
        # type: () -> Any
        attr = '_' + sys._getframe().f_code.co_name
        if not hasattr(self, attr):
            if self.Parser is not CParser:
                setattr(self, attr, self.Parser(loader=self))
            else:
                if getattr(self, '_stream', None) is None:
                    # wait for the stream
                    return None
                else:
                    # if not hasattr(self._stream, 'read') and hasattr(self._stream, 'open'):
                    #     # pathlib.Path() instance
                    #     setattr(self, attr, CParser(self._stream))
                    # else:
                    setattr(self, attr, CParser(self._stream))
                    # self._parser = self._composer = self
                    # nprint('scanner', self.loader.scanner)

        return getattr(self, attr)

    @property
    def composer(self):
        # type: () -> Any
        attr = '_' + sys._getframe().f_code.co_name
        if not hasattr(self, attr):
            setattr(self, attr, self.Composer(loader=self))
        return getattr(self, attr)

    @property
    def constructor(self):
        # type: () -> Any
        attr = '_' + sys._getframe().f_code.co_name
        if not hasattr(self, attr):
            cnst = self.Constructor(preserve_quotes=self.preserve_quotes, loader=self)
            cnst.allow_duplicate_keys = self.allow_duplicate_keys
            setattr(self, attr, cnst)
        return getattr(self, attr)

    @property
    def resolver(self):
        # type: () -> Any
        attr = '_' + sys._getframe().f_code.co_name
        if not hasattr(self, attr):
            setattr(self, attr, self.Resolver(version=self.version, loader=self))
        return getattr(self, attr)

    @property
    def emitter(self):
        # type: () -> Any
        attr = '_' + sys._getframe().f_code.co_name
        if not hasattr(self, attr):
            if self.Emitter is not CEmitter:
                _emitter = self.Emitter(
                    None,
                    canonical=self.canonical,
                    indent=self.old_indent,
                    width=self.width,
                    allow_unicode=self.allow_unicode,
                    line_break=self.line_break,
                    prefix_colon=self.prefix_colon,
                    brace_single_entry_mapping_in_flow_sequence=self.brace_single_entry_mapping_in_flow_sequence,  # NOQA
                    dumper=self,
                )
                setattr(self, attr, _emitter)
                if self.map_indent is not None:
                    _emitter.best_map_indent = self.map_indent
                if self.sequence_indent is not None:
                    _emitter.best_sequence_indent = self.sequence_indent
                if self.sequence_dash_offset is not None:
                    _emitter.sequence_dash_offset = self.sequence_dash_offset
                    # _emitter.block_seq_indent = self.sequence_dash_offset
                if self.compact_seq_seq is not None:
                    _emitter.compact_seq_seq = self.compact_seq_seq
                if self.compact_seq_map is not None:
                    _emitter.compact_seq_map = self.compact_seq_map
            else:
                if getattr(self, '_stream', None) is None:
                    # wait for the stream
                    return None
                return None
        return getattr(self, attr)

    @property
    def serializer(self):
        # type: () -> Any
        attr = '_' + sys._getframe().f_code.co_name
        if not hasattr(self, attr):
            setattr(
                self,
                attr,
                self.Serializer(
                    encoding=self.encoding,
                    explicit_start=self.explicit_start,
                    explicit_end=self.explicit_end,
                    version=self.version,
                    tags=self.tags,
                    dumper=self,
                ),
            )
        return getattr(self, attr)

    @property
    def representer(self):
        # type: () -> Any
        attr = '_' + sys._getframe().f_code.co_name
        if not hasattr(self, attr):
            repres = self.Representer(
                default_style=self.default_style,
                default_flow_style=self.default_flow_style,
                dumper=self,
            )
            if self.sort_base_mapping_type_on_output is not None:
                repres.sort_base_mapping_type_on_output = self.sort_base_mapping_type_on_output
            setattr(self, attr, repres)
        return getattr(self, attr)

    # separate output resolver?

    # def load(self, stream=None):
    #     if self._context_manager:
    #        if not self._input:
    #             raise TypeError("Missing input stream while dumping from context manager")
    #         for data in self._context_manager.load():
    #             yield data
    #         return
    #     if stream is None:
    #         raise TypeError("Need a stream argument when not loading from context manager")
    #     return self.load_one(stream)

    def load(self, stream):
        # type: (Union[Path, StreamTextType]) -> Any
        """
        at this point you either have the non-pure Parser (which has its own reader and
        scanner) or you have the pure Parser.
        If the pure Parser is set, then set the Reader and Scanner, if not already set.
        If either the Scanner or Reader are set, you cannot use the non-pure Parser,
            so reset it to the pure parser and set the Reader resp. Scanner if necessary
        """
        if not hasattr(stream, 'read') and hasattr(stream, 'open'):
            # pathlib.Path() instance
            with stream.open('rb') as fp:  # type: ignore
                return self.load(fp)
        constructor, parser = self.get_constructor_parser(stream)
        try:
            return constructor.get_single_data()
        finally:
            parser.dispose()
            try:
                self._reader.reset_reader()
            except AttributeError:
                pass
            try:
                self._scanner.reset_scanner()
            except AttributeError:
                pass

    def load_all(self, stream, _kw=enforce):  # , skip=None):
        # type: (Union[Path, StreamTextType], Any) -> Any
        if _kw is not enforce:
            raise TypeError(
                '{}.__init__() takes no positional argument but at least '
                'one was given ({!r})'.format(self.__class__.__name__, _kw)
            )
        if not hasattr(stream, 'read') and hasattr(stream, 'open'):
            # pathlib.Path() instance
            with stream.open('r') as fp:  # type: ignore
                for d in self.load_all(fp, _kw=enforce):
                    yield d
                return
        # if skip is None:
        #     skip = []
        # elif isinstance(skip, int):
        #     skip = [skip]
        constructor, parser = self.get_constructor_parser(stream)
        try:
            while constructor.check_data():
                yield constructor.get_data()
        finally:
            parser.dispose()
            try:
                self._reader.reset_reader()
            except AttributeError:
                pass
            try:
                self._scanner.reset_scanner()
            except AttributeError:
                pass

    def get_constructor_parser(self, stream):
        # type: (StreamTextType) -> Any
        """
        the old cyaml needs special setup, and therefore the stream
        """
        if self.Parser is not CParser:
            if self.Reader is None:
                self.Reader = ruamel.yaml.reader.Reader
            if self.Scanner is None:
                self.Scanner = ruamel.yaml.scanner.Scanner
            self.reader.stream = stream
        else:
            if self.Reader is not None:
                if self.Scanner is None:
                    self.Scanner = ruamel.yaml.scanner.Scanner
                self.Parser = ruamel.yaml.parser.Parser
                self.reader.stream = stream
            elif self.Scanner is not None:
                if self.Reader is None:
                    self.Reader = ruamel.yaml.reader.Reader
                self.Parser = ruamel.yaml.parser.Parser
                self.reader.stream = stream
            else:
                # combined C level reader>scanner>parser
                # does some calls to the resolver, e.g. BaseResolver.descend_resolver
                # if you just initialise the CParser, to much of resolver.py
                # is actually used
                rslvr = self.Resolver
                # if rslvr is ruamel.yaml.resolver.VersionedResolver:
                #     rslvr = ruamel.yaml.resolver.Resolver

                class XLoader(self.Parser, self.Constructor, rslvr):  # type: ignore
                    def __init__(selfx, stream, version=self.version, preserve_quotes=None):
                        # type: (StreamTextType, Optional[VersionType], Optional[bool]) -> None  # NOQA
                        CParser.__init__(selfx, stream)
                        selfx._parser = selfx._composer = selfx
                        self.Constructor.__init__(selfx, loader=selfx)
                        selfx.allow_duplicate_keys = self.allow_duplicate_keys
                        rslvr.__init__(selfx, version=version, loadumper=selfx)

                self._stream = stream
                loader = XLoader(stream)
                return loader, loader
        return self.constructor, self.parser

    def dump(self, data, stream=None, _kw=enforce, transform=None):
        # type: (Any, Union[Path, StreamType], Any, Any) -> Any
        if self._context_manager:
            if not self._output:
                raise TypeError('Missing output stream while dumping from context manager')
            if _kw is not enforce:
                raise TypeError(
                    '{}.dump() takes one positional argument but at least '
                    'two were given ({!r})'.format(self.__class__.__name__, _kw)
                )
            if transform is not None:
                raise TypeError(
                    '{}.dump() in the context manager cannot have transform keyword '
                    ''.format(self.__class__.__name__)
                )
            self._context_manager.dump(data)
        else:  # old style
            if stream is None:
                raise TypeError('Need a stream argument when not dumping from context manager')
            return self.dump_all([data], stream, _kw, transform=transform)

    def dump_all(self, documents, stream, _kw=enforce, transform=None):
        # type: (Any, Union[Path, StreamType], Any, Any) -> Any
        if self._context_manager:
            raise NotImplementedError
        if _kw is not enforce:
            raise TypeError(
                '{}.dump(_all) takes two positional argument but at least '
                'three were given ({!r})'.format(self.__class__.__name__, _kw)
            )
        self._output = stream
        self._context_manager = YAMLContextManager(self, transform=transform)
        for data in documents:
            self._context_manager.dump(data)
        self._context_manager.teardown_output()
        self._output = None
        self._context_manager = None

    def Xdump_all(self, documents, stream, _kw=enforce, transform=None):
        # type: (Any, Union[Path, StreamType], Any, Any) -> Any
        """
        Serialize a sequence of Python objects into a YAML stream.
        """
        if not hasattr(stream, 'write') and hasattr(stream, 'open'):
            # pathlib.Path() instance
            with stream.open('w') as fp:
                return self.dump_all(documents, fp, _kw, transform=transform)
        if _kw is not enforce:
            raise TypeError(
                '{}.dump(_all) takes two positional argument but at least '
                'three were given ({!r})'.format(self.__class__.__name__, _kw)
            )
        # The stream should have the methods `write` and possibly `flush`.
        if self.top_level_colon_align is True:
            tlca = max([len(str(x)) for x in documents[0]])  # type: Any
        else:
            tlca = self.top_level_colon_align
        if transform is not None:
            fstream = stream
            if self.encoding is None:
                stream = StringIO()
            else:
                stream = BytesIO()
        serializer, representer, emitter = self.get_serializer_representer_emitter(
            stream, tlca
        )
        try:
            self.serializer.open()
            for data in documents:
                try:
                    self.representer.represent(data)
                except AttributeError:
                    # nprint(dir(dumper._representer))
                    raise
            self.serializer.close()
        finally:
            try:
                self.emitter.dispose()
            except AttributeError:
                raise
                # self.dumper.dispose()  # cyaml
            delattr(self, '_serializer')
            delattr(self, '_emitter')
        if transform:
            val = stream.getvalue()
            if self.encoding:
                val = val.decode(self.encoding)
            if fstream is None:
                transform(val)
            else:
                fstream.write(transform(val))
        return None

    def get_serializer_representer_emitter(self, stream, tlca):
        # type: (StreamType, Any) -> Any
        # we have only .Serializer to deal with (vs .Reader & .Scanner), much simpler
        if self.Emitter is not CEmitter:
            if self.Serializer is None:
                self.Serializer = ruamel.yaml.serializer.Serializer
            self.emitter.stream = stream
            self.emitter.top_level_colon_align = tlca
            return self.serializer, self.representer, self.emitter
        if self.Serializer is not None:
            # cannot set serializer with CEmitter
            self.Emitter = ruamel.yaml.emitter.Emitter
            self.emitter.stream = stream
            self.emitter.top_level_colon_align = tlca
            return self.serializer, self.representer, self.emitter
        # C routines

        rslvr = (
            ruamel.yaml.resolver.BaseResolver
            if self.typ == 'base'
            else ruamel.yaml.resolver.Resolver
        )

        class XDumper(CEmitter, self.Representer, rslvr):  # type: ignore
            def __init__(
                selfx,
                stream,
                default_style=None,
                default_flow_style=None,
                canonical=None,
                indent=None,
                width=None,
                allow_unicode=None,
                line_break=None,
                encoding=None,
                explicit_start=None,
                explicit_end=None,
                version=None,
                tags=None,
                block_seq_indent=None,
                top_level_colon_align=None,
                prefix_colon=None,
            ):
                # type: (StreamType, Any, Any, Any, Optional[bool], Optional[int], Optional[int], Optional[bool], Any, Any, Optional[bool], Optional[bool], Any, Any, Any, Any, Any) -> None   # NOQA
                CEmitter.__init__(
                    selfx,
                    stream,
                    canonical=canonical,
                    indent=indent,
                    width=width,
                    encoding=encoding,
                    allow_unicode=allow_unicode,
                    line_break=line_break,
                    explicit_start=explicit_start,
                    explicit_end=explicit_end,
                    version=version,
                    tags=tags,
                )
                selfx._emitter = selfx._serializer = selfx._representer = selfx
                self.Representer.__init__(
                    selfx, default_style=default_style, default_flow_style=default_flow_style
                )
                rslvr.__init__(selfx)

        self._stream = stream
        dumper = XDumper(
            stream,
            default_style=self.default_style,
            default_flow_style=self.default_flow_style,
            canonical=self.canonical,
            indent=self.old_indent,
            width=self.width,
            allow_unicode=self.allow_unicode,
            line_break=self.line_break,
            explicit_start=self.explicit_start,
            explicit_end=self.explicit_end,
            version=self.version,
            tags=self.tags,
        )
        self._emitter = self._serializer = dumper
        return dumper, dumper, dumper

    # basic types
    def map(self, **kw):
        # type: (Any) -> Any
        if self.typ == 'rt':
            from ruamel.yaml.comments import CommentedMap

            return CommentedMap(**kw)
        else:
            return dict(**kw)

    def seq(self, *args):
        # type: (Any) -> Any
        if self.typ == 'rt':
            from ruamel.yaml.comments import CommentedSeq

            return CommentedSeq(*args)
        else:
            return list(*args)

    # helpers
    def official_plug_ins(self):
        # type: () -> Any
        bd = os.path.dirname(__file__)
        gpbd = os.path.dirname(os.path.dirname(bd))
        res = [x.replace(gpbd, "")[1:-3] for x in glob.glob(bd + '/*/__plug_in__.py')]
        return res

    def register_class(self, cls):
        # type:(Any) -> Any
        """
        register a class for dumping loading
        - if it has attribute yaml_tag use that to register, else use class name
        - if it has methods to_yaml/from_yaml use those to dump/load else dump attributes
          as mapping
        """
        tag = getattr(cls, 'yaml_tag', '!' + cls.__name__)
        try:
            self.representer.add_representer(cls, cls.to_yaml)
        except AttributeError:

            def t_y(representer, data):
                # type: (Any, Any) -> Any
                return representer.represent_yaml_object(
                    tag, data, cls, flow_style=representer.default_flow_style
                )

            self.representer.add_representer(cls, t_y)
        try:
            self.constructor.add_constructor(tag, cls.from_yaml)
        except AttributeError:

            def f_y(constructor, node):
                # type: (Any, Any) -> Any
                return constructor.construct_yaml_object(node, cls)

            self.constructor.add_constructor(tag, f_y)
        return cls

    def parse(self, stream):
        # type: (StreamTextType) -> Any
        """
        Parse a YAML stream and produce parsing events.
        """
        _, parser = self.get_constructor_parser(stream)
        try:
            while parser.check_event():
                yield parser.get_event()
        finally:
            parser.dispose()
            try:
                self._reader.reset_reader()
            except AttributeError:
                pass
            try:
                self._scanner.reset_scanner()
            except AttributeError:
                pass

    # ### context manager

    def __enter__(self):
        # type: () -> Any
        self._context_manager = YAMLContextManager(self)
        return self

    def __exit__(self, typ, value, traceback):
        # type: (Any, Any, Any) -> None
        if typ:
            nprint('typ', typ)
        self._context_manager.teardown_output()
        # self._context_manager.teardown_input()
        self._context_manager = None

    # ### backwards compatibility
    def _indent(self, mapping=None, sequence=None, offset=None):
        # type: (Any, Any, Any) -> None
        if mapping is not None:
            self.map_indent = mapping
        if sequence is not None:
            self.sequence_indent = sequence
        if offset is not None:
            self.sequence_dash_offset = offset

    @property
    def indent(self):
        # type: () -> Any
        return self._indent

    @indent.setter
    def indent(self, val):
        # type: (Any) -> None
        self.old_indent = val

    @property
    def block_seq_indent(self):
        # type: () -> Any
        return self.sequence_dash_offset

    @block_seq_indent.setter
    def block_seq_indent(self, val):
        # type: (Any) -> None
        self.sequence_dash_offset = val

    def compact(self, seq_seq=None, seq_map=None):
        # type: (Any, Any) -> None
        self.compact_seq_seq = seq_seq
        self.compact_seq_map = seq_map


class YAMLContextManager(object):
    def __init__(self, yaml, transform=None):
        # type: (Any, Any) -> None  # used to be: (Any, Optional[Callable]) -> None
        self._yaml = yaml
        self._output_inited = False
        self._output_path = None
        self._output = self._yaml._output
        self._transform = transform

        # self._input_inited = False
        # self._input = input
        # self._input_path = None
        # self._transform = yaml.transform
        # self._fstream = None

        if not hasattr(self._output, 'write') and hasattr(self._output, 'open'):
            # pathlib.Path() instance, open with the same mode
            self._output_path = self._output
            self._output = self._output_path.open('w')

        # if not hasattr(self._stream, 'write') and hasattr(stream, 'open'):
        # if not hasattr(self._input, 'read') and hasattr(self._input, 'open'):
        #    # pathlib.Path() instance, open with the same mode
        #    self._input_path = self._input
        #    self._input = self._input_path.open('r')

        if self._transform is not None:
            self._fstream = self._output
            if self._yaml.encoding is None:
                self._output = StringIO()
            else:
                self._output = BytesIO()

    def teardown_output(self):
        # type: () -> None
        if self._output_inited:
            self._yaml.serializer.close()
        else:
            return
        try:
            self._yaml.emitter.dispose()
        except AttributeError:
            raise
            # self.dumper.dispose()  # cyaml
        try:
            delattr(self._yaml, '_serializer')
            delattr(self._yaml, '_emitter')
        except AttributeError:
            raise
        if self._transform:
            val = self._output.getvalue()
            if self._yaml.encoding:
                val = val.decode(self._yaml.encoding)
            if self._fstream is None:
                self._transform(val)
            else:
                self._fstream.write(self._transform(val))
                self._fstream.flush()
                self._output = self._fstream  # maybe not necessary
        if self._output_path is not None:
            self._output.close()

    def init_output(self, first_data):
        # type: (Any) -> None
        if self._yaml.top_level_colon_align is True:
            tlca = max([len(str(x)) for x in first_data])  # type: Any
        else:
            tlca = self._yaml.top_level_colon_align
        self._yaml.get_serializer_representer_emitter(self._output, tlca)
        self._yaml.serializer.open()
        self._output_inited = True

    def dump(self, data):
        # type: (Any) -> None
        if not self._output_inited:
            self.init_output(data)
        try:
            self._yaml.representer.represent(data)
        except AttributeError:
            # nprint(dir(dumper._representer))
            raise

    # def teardown_input(self):
    #     pass
    #
    # def init_input(self):
    #     # set the constructor and parser on YAML() instance
    #     self._yaml.get_constructor_parser(stream)
    #
    # def load(self):
    #     if not self._input_inited:
    #         self.init_input()
    #     try:
    #         while self._yaml.constructor.check_data():
    #             yield self._yaml.constructor.get_data()
    #     finally:
    #         parser.dispose()
    #         try:
    #             self._reader.reset_reader()  # type: ignore
    #         except AttributeError:
    #             pass
    #         try:
    #             self._scanner.reset_scanner()  # type: ignore
    #         except AttributeError:
    #             pass


def yaml_object(yml):
    # type: (Any) -> Any
    """ decorator for classes that needs to dump/load objects
    The tag for such objects is taken from the class attribute yaml_tag (or the
    class name in lowercase in case unavailable)
    If methods to_yaml and/or from_yaml are available, these are called for dumping resp.
    loading, default routines (dumping a mapping of the attributes) used otherwise.
    """

    def yo_deco(cls):
        # type: (Any) -> Any
        tag = getattr(cls, 'yaml_tag', '!' + cls.__name__)
        try:
            yml.representer.add_representer(cls, cls.to_yaml)
        except AttributeError:

            def t_y(representer, data):
                # type: (Any, Any) -> Any
                return representer.represent_yaml_object(
                    tag, data, cls, flow_style=representer.default_flow_style
                )

            yml.representer.add_representer(cls, t_y)
        try:
            yml.constructor.add_constructor(tag, cls.from_yaml)
        except AttributeError:

            def f_y(constructor, node):
                # type: (Any, Any) -> Any
                return constructor.construct_yaml_object(node, cls)

            yml.constructor.add_constructor(tag, f_y)
        return cls

    return yo_deco


########################################################################################


def scan(stream, Loader=Loader):
    # type: (StreamTextType, Any) -> Any
    """
    Scan a YAML stream and produce scanning tokens.
    """
    loader = Loader(stream)
    try:
        while loader.scanner.check_token():
            yield loader.scanner.get_token()
    finally:
        loader._parser.dispose()


def parse(stream, Loader=Loader):
    # type: (StreamTextType, Any) -> Any
    """
    Parse a YAML stream and produce parsing events.
    """
    loader = Loader(stream)
    try:
        while loader._parser.check_event():
            yield loader._parser.get_event()
    finally:
        loader._parser.dispose()


def compose(stream, Loader=Loader):
    # type: (StreamTextType, Any) -> Any
    """
    Parse the first YAML document in a stream
    and produce the corresponding representation tree.
    """
    loader = Loader(stream)
    try:
        return loader.get_single_node()
    finally:
        loader.dispose()


def compose_all(stream, Loader=Loader):
    # type: (StreamTextType, Any) -> Any
    """
    Parse all YAML documents in a stream
    and produce corresponding representation trees.
    """
    loader = Loader(stream)
    try:
        while loader.check_node():
            yield loader._composer.get_node()
    finally:
        loader._parser.dispose()


def load(stream, Loader=None, version=None, preserve_quotes=None):
    # type: (StreamTextType, Any, Optional[VersionType], Any) -> Any
    """
    Parse the first YAML document in a stream
    and produce the corresponding Python object.
    """
    if Loader is None:
        warnings.warn(UnsafeLoaderWarning.text, UnsafeLoaderWarning, stacklevel=2)
        Loader = UnsafeLoader
    loader = Loader(stream, version, preserve_quotes=preserve_quotes)
    try:
        return loader._constructor.get_single_data()
    finally:
        loader._parser.dispose()
        try:
            loader._reader.reset_reader()
        except AttributeError:
            pass
        try:
            loader._scanner.reset_scanner()
        except AttributeError:
            pass


def load_all(stream, Loader=None, version=None, preserve_quotes=None):
    # type: (Optional[StreamTextType], Any, Optional[VersionType], Optional[bool]) -> Any  # NOQA
    """
    Parse all YAML documents in a stream
    and produce corresponding Python objects.
    """
    if Loader is None:
        warnings.warn(UnsafeLoaderWarning.text, UnsafeLoaderWarning, stacklevel=2)
        Loader = UnsafeLoader
    loader = Loader(stream, version, preserve_quotes=preserve_quotes)
    try:
        while loader._constructor.check_data():
            yield loader._constructor.get_data()
    finally:
        loader._parser.dispose()
        try:
            loader._reader.reset_reader()
        except AttributeError:
            pass
        try:
            loader._scanner.reset_scanner()
        except AttributeError:
            pass


def safe_load(stream, version=None):
    # type: (StreamTextType, Optional[VersionType]) -> Any
    """
    Parse the first YAML document in a stream
    and produce the corresponding Python object.
    Resolve only basic YAML tags.
    """
    return load(stream, SafeLoader, version)


def safe_load_all(stream, version=None):
    # type: (StreamTextType, Optional[VersionType]) -> Any
    """
    Parse all YAML documents in a stream
    and produce corresponding Python objects.
    Resolve only basic YAML tags.
    """
    return load_all(stream, SafeLoader, version)


def round_trip_load(stream, version=None, preserve_quotes=None):
    # type: (StreamTextType, Optional[VersionType], Optional[bool]) -> Any
    """
    Parse the first YAML document in a stream
    and produce the corresponding Python object.
    Resolve only basic YAML tags.
    """
    return load(stream, RoundTripLoader, version, preserve_quotes=preserve_quotes)


def round_trip_load_all(stream, version=None, preserve_quotes=None):
    # type: (StreamTextType, Optional[VersionType], Optional[bool]) -> Any
    """
    Parse all YAML documents in a stream
    and produce corresponding Python objects.
    Resolve only basic YAML tags.
    """
    return load_all(stream, RoundTripLoader, version, preserve_quotes=preserve_quotes)


def emit(
    events,
    stream=None,
    Dumper=Dumper,
    canonical=None,
    indent=None,
    width=None,
    allow_unicode=None,
    line_break=None,
):
    # type: (Any, Optional[StreamType], Any, Optional[bool], Union[int, None], Optional[int], Optional[bool], Any) -> Any  # NOQA
    """
    Emit YAML parsing events into a stream.
    If stream is None, return the produced string instead.
    """
    getvalue = None
    if stream is None:
        stream = StringIO()
        getvalue = stream.getvalue
    dumper = Dumper(
        stream,
        canonical=canonical,
        indent=indent,
        width=width,
        allow_unicode=allow_unicode,
        line_break=line_break,
    )
    try:
        for event in events:
            dumper.emit(event)
    finally:
        try:
            dumper._emitter.dispose()
        except AttributeError:
            raise
            dumper.dispose()  # cyaml
    if getvalue is not None:
        return getvalue()


enc = None if PY3 else 'utf-8'


def serialize_all(
    nodes,
    stream=None,
    Dumper=Dumper,
    canonical=None,
    indent=None,
    width=None,
    allow_unicode=None,
    line_break=None,
    encoding=enc,
    explicit_start=None,
    explicit_end=None,
    version=None,
    tags=None,
):
    # type: (Any, Optional[StreamType], Any, Any, Optional[int], Optional[int], Optional[bool], Any, Any, Optional[bool], Optional[bool], Optional[VersionType], Any) -> Any # NOQA
    """
    Serialize a sequence of representation trees into a YAML stream.
    If stream is None, return the produced string instead.
    """
    getvalue = None
    if stream is None:
        if encoding is None:
            stream = StringIO()
        else:
            stream = BytesIO()
        getvalue = stream.getvalue
    dumper = Dumper(
        stream,
        canonical=canonical,
        indent=indent,
        width=width,
        allow_unicode=allow_unicode,
        line_break=line_break,
        encoding=encoding,
        version=version,
        tags=tags,
        explicit_start=explicit_start,
        explicit_end=explicit_end,
    )
    try:
        dumper._serializer.open()
        for node in nodes:
            dumper.serialize(node)
        dumper._serializer.close()
    finally:
        try:
            dumper._emitter.dispose()
        except AttributeError:
            raise
            dumper.dispose()  # cyaml
    if getvalue is not None:
        return getvalue()


def serialize(node, stream=None, Dumper=Dumper, **kwds):
    # type: (Any, Optional[StreamType], Any, Any) -> Any
    """
    Serialize a representation tree into a YAML stream.
    If stream is None, return the produced string instead.
    """
    return serialize_all([node], stream, Dumper=Dumper, **kwds)


def dump_all(
    documents,
    stream=None,
    Dumper=Dumper,
    default_style=None,
    default_flow_style=None,
    canonical=None,
    indent=None,
    width=None,
    allow_unicode=None,
    line_break=None,
    encoding=enc,
    explicit_start=None,
    explicit_end=None,
    version=None,
    tags=None,
    block_seq_indent=None,
    top_level_colon_align=None,
    prefix_colon=None,
):
    # type: (Any, Optional[StreamType], Any, Any, Any, Optional[bool], Optional[int], Optional[int], Optional[bool], Any, Any, Optional[bool], Optional[bool], Any, Any, Any, Any, Any) -> Optional[str]   # NOQA
    """
    Serialize a sequence of Python objects into a YAML stream.
    If stream is None, return the produced string instead.
    """
    getvalue = None
    if top_level_colon_align is True:
        top_level_colon_align = max([len(str(x)) for x in documents[0]])
    if stream is None:
        if encoding is None:
            stream = StringIO()
        else:
            stream = BytesIO()
        getvalue = stream.getvalue
    dumper = Dumper(
        stream,
        default_style=default_style,
        default_flow_style=default_flow_style,
        canonical=canonical,
        indent=indent,
        width=width,
        allow_unicode=allow_unicode,
        line_break=line_break,
        encoding=encoding,
        explicit_start=explicit_start,
        explicit_end=explicit_end,
        version=version,
        tags=tags,
        block_seq_indent=block_seq_indent,
        top_level_colon_align=top_level_colon_align,
        prefix_colon=prefix_colon,
    )
    try:
        dumper._serializer.open()
        for data in documents:
            try:
                dumper._representer.represent(data)
            except AttributeError:
                # nprint(dir(dumper._representer))
                raise
        dumper._serializer.close()
    finally:
        try:
            dumper._emitter.dispose()
        except AttributeError:
            raise
            dumper.dispose()  # cyaml
    if getvalue is not None:
        return getvalue()
    return None


def dump(
    data,
    stream=None,
    Dumper=Dumper,
    default_style=None,
    default_flow_style=None,
    canonical=None,
    indent=None,
    width=None,
    allow_unicode=None,
    line_break=None,
    encoding=enc,
    explicit_start=None,
    explicit_end=None,
    version=None,
    tags=None,
    block_seq_indent=None,
):
    # type: (Any, Optional[StreamType], Any, Any, Any, Optional[bool], Optional[int], Optional[int], Optional[bool], Any, Any, Optional[bool], Optional[bool], Optional[VersionType], Any, Any) -> Optional[str]   # NOQA
    """
    Serialize a Python object into a YAML stream.
    If stream is None, return the produced string instead.

    default_style ∈ None, '', '"', "'", '|', '>'

    """
    return dump_all(
        [data],
        stream,
        Dumper=Dumper,
        default_style=default_style,
        default_flow_style=default_flow_style,
        canonical=canonical,
        indent=indent,
        width=width,
        allow_unicode=allow_unicode,
        line_break=line_break,
        encoding=encoding,
        explicit_start=explicit_start,
        explicit_end=explicit_end,
        version=version,
        tags=tags,
        block_seq_indent=block_seq_indent,
    )


def safe_dump_all(documents, stream=None, **kwds):
    # type: (Any, Optional[StreamType], Any) -> Optional[str]
    """
    Serialize a sequence of Python objects into a YAML stream.
    Produce only basic YAML tags.
    If stream is None, return the produced string instead.
    """
    return dump_all(documents, stream, Dumper=SafeDumper, **kwds)


def safe_dump(data, stream=None, **kwds):
    # type: (Any, Optional[StreamType], Any) -> Optional[str]
    """
    Serialize a Python object into a YAML stream.
    Produce only basic YAML tags.
    If stream is None, return the produced string instead.
    """
    return dump_all([data], stream, Dumper=SafeDumper, **kwds)


def round_trip_dump(
    data,
    stream=None,
    Dumper=RoundTripDumper,
    default_style=None,
    default_flow_style=None,
    canonical=None,
    indent=None,
    width=None,
    allow_unicode=None,
    line_break=None,
    encoding=enc,
    explicit_start=None,
    explicit_end=None,
    version=None,
    tags=None,
    block_seq_indent=None,
    top_level_colon_align=None,
    prefix_colon=None,
):
    # type: (Any, Optional[StreamType], Any, Any, Any, Optional[bool], Optional[int], Optional[int], Optional[bool], Any, Any, Optional[bool], Optional[bool], Optional[VersionType], Any, Any, Any, Any) -> Optional[str]   # NOQA
    allow_unicode = True if allow_unicode is None else allow_unicode
    return dump_all(
        [data],
        stream,
        Dumper=Dumper,
        default_style=default_style,
        default_flow_style=default_flow_style,
        canonical=canonical,
        indent=indent,
        width=width,
        allow_unicode=allow_unicode,
        line_break=line_break,
        encoding=encoding,
        explicit_start=explicit_start,
        explicit_end=explicit_end,
        version=version,
        tags=tags,
        block_seq_indent=block_seq_indent,
        top_level_colon_align=top_level_colon_align,
        prefix_colon=prefix_colon,
    )


# Loader/Dumper are no longer composites, to get to the associated
# Resolver()/Representer(), etc., you need to instantiate the class


def add_implicit_resolver(
    tag, regexp, first=None, Loader=None, Dumper=None, resolver=Resolver
):
    # type: (Any, Any, Any, Any, Any, Any) -> None
    """
    Add an implicit scalar detector.
    If an implicit scalar value matches the given regexp,
    the corresponding tag is assigned to the scalar.
    first is a sequence of possible initial characters or None.
    """
    if Loader is None and Dumper is None:
        resolver.add_implicit_resolver(tag, regexp, first)
        return
    if Loader:
        if hasattr(Loader, 'add_implicit_resolver'):
            Loader.add_implicit_resolver(tag, regexp, first)
        elif issubclass(
            Loader, (BaseLoader, SafeLoader, ruamel.yaml.loader.Loader, RoundTripLoader)
        ):
            Resolver.add_implicit_resolver(tag, regexp, first)
        else:
            raise NotImplementedError
    if Dumper:
        if hasattr(Dumper, 'add_implicit_resolver'):
            Dumper.add_implicit_resolver(tag, regexp, first)
        elif issubclass(
            Dumper, (BaseDumper, SafeDumper, ruamel.yaml.dumper.Dumper, RoundTripDumper)
        ):
            Resolver.add_implicit_resolver(tag, regexp, first)
        else:
            raise NotImplementedError


# this code currently not tested
def add_path_resolver(tag, path, kind=None, Loader=None, Dumper=None, resolver=Resolver):
    # type: (Any, Any, Any, Any, Any, Any) -> None
    """
    Add a path based resolver for the given tag.
    A path is a list of keys that forms a path
    to a node in the representation tree.
    Keys can be string values, integers, or None.
    """
    if Loader is None and Dumper is None:
        resolver.add_path_resolver(tag, path, kind)
        return
    if Loader:
        if hasattr(Loader, 'add_path_resolver'):
            Loader.add_path_resolver(tag, path, kind)
        elif issubclass(
            Loader, (BaseLoader, SafeLoader, ruamel.yaml.loader.Loader, RoundTripLoader)
        ):
            Resolver.add_path_resolver(tag, path, kind)
        else:
            raise NotImplementedError
    if Dumper:
        if hasattr(Dumper, 'add_path_resolver'):
            Dumper.add_path_resolver(tag, path, kind)
        elif issubclass(
            Dumper, (BaseDumper, SafeDumper, ruamel.yaml.dumper.Dumper, RoundTripDumper)
        ):
            Resolver.add_path_resolver(tag, path, kind)
        else:
            raise NotImplementedError


def add_constructor(tag, object_constructor, Loader=None, constructor=Constructor):
    # type: (Any, Any, Any, Any) -> None
    """
    Add an object constructor for the given tag.
    object_onstructor is a function that accepts a Loader instance
    and a node object and produces the corresponding Python object.
    """
    if Loader is None:
        constructor.add_constructor(tag, object_constructor)
    else:
        if hasattr(Loader, 'add_constructor'):
            Loader.add_constructor(tag, object_constructor)
            return
        if issubclass(Loader, BaseLoader):
            BaseConstructor.add_constructor(tag, object_constructor)
        elif issubclass(Loader, SafeLoader):
            SafeConstructor.add_constructor(tag, object_constructor)
        elif issubclass(Loader, Loader):
            Constructor.add_constructor(tag, object_constructor)
        elif issubclass(Loader, RoundTripLoader):
            RoundTripConstructor.add_constructor(tag, object_constructor)
        else:
            raise NotImplementedError


def add_multi_constructor(tag_prefix, multi_constructor, Loader=None, constructor=Constructor):
    # type: (Any, Any, Any, Any) -> None
    """
    Add a multi-constructor for the given tag prefix.
    Multi-constructor is called for a node if its tag starts with tag_prefix.
    Multi-constructor accepts a Loader instance, a tag suffix,
    and a node object and produces the corresponding Python object.
    """
    if Loader is None:
        constructor.add_multi_constructor(tag_prefix, multi_constructor)
    else:
        if False and hasattr(Loader, 'add_multi_constructor'):
            Loader.add_multi_constructor(tag_prefix, constructor)
            return
        if issubclass(Loader, BaseLoader):
            BaseConstructor.add_multi_constructor(tag_prefix, multi_constructor)
        elif issubclass(Loader, SafeLoader):
            SafeConstructor.add_multi_constructor(tag_prefix, multi_constructor)
        elif issubclass(Loader, ruamel.yaml.loader.Loader):
            Constructor.add_multi_constructor(tag_prefix, multi_constructor)
        elif issubclass(Loader, RoundTripLoader):
            RoundTripConstructor.add_multi_constructor(tag_prefix, multi_constructor)
        else:
            raise NotImplementedError


def add_representer(data_type, object_representer, Dumper=None, representer=Representer):
    # type: (Any, Any, Any, Any) -> None
    """
    Add a representer for the given type.
    object_representer is a function accepting a Dumper instance
    and an instance of the given data type
    and producing the corresponding representation node.
    """
    if Dumper is None:
        representer.add_representer(data_type, object_representer)
    else:
        if hasattr(Dumper, 'add_representer'):
            Dumper.add_representer(data_type, object_representer)
            return
        if issubclass(Dumper, BaseDumper):
            BaseRepresenter.add_representer(data_type, object_representer)
        elif issubclass(Dumper, SafeDumper):
            SafeRepresenter.add_representer(data_type, object_representer)
        elif issubclass(Dumper, Dumper):
            Representer.add_representer(data_type, object_representer)
        elif issubclass(Dumper, RoundTripDumper):
            RoundTripRepresenter.add_representer(data_type, object_representer)
        else:
            raise NotImplementedError


# this code currently not tested
def add_multi_representer(data_type, multi_representer, Dumper=None, representer=Representer):
    # type: (Any, Any, Any, Any) -> None
    """
    Add a representer for the given type.
    multi_representer is a function accepting a Dumper instance
    and an instance of the given data type or subtype
    and producing the corresponding representation node.
    """
    if Dumper is None:
        representer.add_multi_representer(data_type, multi_representer)
    else:
        if hasattr(Dumper, 'add_multi_representer'):
            Dumper.add_multi_representer(data_type, multi_representer)
            return
        if issubclass(Dumper, BaseDumper):
            BaseRepresenter.add_multi_representer(data_type, multi_representer)
        elif issubclass(Dumper, SafeDumper):
            SafeRepresenter.add_multi_representer(data_type, multi_representer)
        elif issubclass(Dumper, Dumper):
            Representer.add_multi_representer(data_type, multi_representer)
        elif issubclass(Dumper, RoundTripDumper):
            RoundTripRepresenter.add_multi_representer(data_type, multi_representer)
        else:
            raise NotImplementedError


class YAMLObjectMetaclass(type):
    """
    The metaclass for YAMLObject.
    """

    def __init__(cls, name, bases, kwds):
        # type: (Any, Any, Any) -> None
        super(YAMLObjectMetaclass, cls).__init__(name, bases, kwds)
        if 'yaml_tag' in kwds and kwds['yaml_tag'] is not None:
            cls.yaml_constructor.add_constructor(cls.yaml_tag, cls.from_yaml)  # type: ignore
            cls.yaml_representer.add_representer(cls, cls.to_yaml)  # type: ignore


class YAMLObject(with_metaclass(YAMLObjectMetaclass)):  # type: ignore
    """
    An object that can dump itself to a YAML stream
    and load itself from a YAML stream.
    """

    __slots__ = ()  # no direct instantiation, so allow immutable subclasses

    yaml_constructor = Constructor
    yaml_representer = Representer

    yaml_tag = None  # type: Any
    yaml_flow_style = None  # type: Any

    @classmethod
    def from_yaml(cls, constructor, node):
        # type: (Any, Any) -> Any
        """
        Convert a representation node to a Python object.
        """
        return constructor.construct_yaml_object(node, cls)

    @classmethod
    def to_yaml(cls, representer, data):
        # type: (Any, Any) -> Any
        """
        Convert a Python object to a representation node.
        """
        return representer.represent_yaml_object(
            cls.yaml_tag, data, cls, flow_style=cls.yaml_flow_style
        )
