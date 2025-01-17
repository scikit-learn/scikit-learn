from __future__ import annotations

from typing import TYPE_CHECKING, Any

from sphinx.domains.c._ast import (
    ASTDeclaration,
    ASTIdentifier,
    ASTNestedName,
)
from sphinx.locale import __
from sphinx.util import logging

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Iterator, Sequence

    from typing_extensions import Self

    from sphinx.environment import BuildEnvironment

logger = logging.getLogger(__name__)


class _DuplicateSymbolError(Exception):
    def __init__(self, symbol: Symbol, declaration: ASTDeclaration) -> None:
        assert symbol
        assert declaration
        self.symbol = symbol
        self.declaration = declaration

    def __str__(self) -> str:
        return "Internal C duplicate symbol error:\n%s" % self.symbol.dump(0)


class SymbolLookupResult:
    def __init__(self, symbols: Sequence[Symbol], parentSymbol: Symbol,
                 ident: ASTIdentifier) -> None:
        self.symbols = symbols
        self.parentSymbol = parentSymbol
        self.ident = ident


class LookupKey:
    def __init__(self, data: list[tuple[ASTIdentifier, str]]) -> None:
        self.data = data

    def __str__(self) -> str:
        inner = ', '.join(f"({ident}, {id_})" for ident, id_ in self.data)
        return f'[{inner}]'


class Symbol:
    debug_indent = 0
    debug_indent_string = "  "
    debug_lookup = False
    debug_show_tree = False

    def __copy__(self) -> Self:
        raise AssertionError  # shouldn't happen

    def __deepcopy__(self, memo: Any) -> Symbol:
        if self.parent:
            raise AssertionError  # shouldn't happen
        # the domain base class makes a copy of the initial data, which is fine
        return Symbol(None, None, None, None, None)

    @staticmethod
    def debug_print(*args: Any) -> None:
        msg = Symbol.debug_indent_string * Symbol.debug_indent
        msg += "".join(str(e) for e in args)
        logger.debug(msg)

    def _assert_invariants(self) -> None:
        if not self.parent:
            # parent == None means global scope, so declaration means a parent
            assert not self.declaration
            assert not self.docname
        else:
            if self.declaration:
                assert self.docname

    def __setattr__(self, key: str, value: Any) -> None:
        if key == "children":
            raise AssertionError
        return super().__setattr__(key, value)

    def __init__(
        self,
        parent: Symbol | None,
        ident: ASTIdentifier | None,
        declaration: ASTDeclaration | None,
        docname: str | None,
        line: int | None,
    ) -> None:
        self.parent = parent
        # declarations in a single directive are linked together
        self.siblingAbove: Symbol | None = None
        self.siblingBelow: Symbol | None = None
        self.ident = ident
        self.declaration = declaration
        self.docname = docname
        self.line = line
        self.isRedeclaration = False
        self._assert_invariants()

        # These properties store the same children for different access patterns.
        # ``_add_child()`` and ``_remove_child()`` should be used for modifying them.
        self._children_by_name: dict[str, Symbol] = {}
        self._children_by_docname: dict[str, dict[str, Symbol]] = {}
        self._anon_children: set[Symbol] = set()

        if self.parent:
            self.parent._add_child(self)
        if self.declaration:
            self.declaration.symbol = self

        # Do symbol addition after self._children has been initialised.
        self._add_function_params()

    def __repr__(self) -> str:
        return f'<Symbol {self.to_string(indent=0)!r}>'

    @property
    def _children(self) -> Iterable[Symbol]:
        return self._children_by_name.values()

    def _add_child(self, child: Symbol) -> None:
        name = child.ident.name
        if name in self._children_by_name:
            # Duplicate so don't add - will be reported in _add_symbols()
            return
        self._children_by_name[name] = child
        self._children_by_docname.setdefault(child.docname, {})[name] = child
        if child.ident.is_anonymous:
            self._anon_children.add(child)

    def _remove_child(self, child: Symbol) -> None:
        name = child.ident.name
        self._children_by_name.pop(name, None)
        self._children_by_docname.get(child.docname, {}).pop(name, None)
        if child.ident.is_anonymous:
            self._anon_children.discard(child)

    def _fill_empty(self, declaration: ASTDeclaration, docname: str, line: int) -> None:
        self._assert_invariants()
        assert self.declaration is None
        assert self.docname is None
        assert self.line is None
        assert declaration is not None
        assert docname is not None
        assert line is not None
        self.declaration = declaration
        self.declaration.symbol = self
        self.docname = docname
        self.line = line
        self._assert_invariants()
        # and symbol addition should be done as well
        self._add_function_params()

    def _add_function_params(self) -> None:
        if Symbol.debug_lookup:
            Symbol.debug_indent += 1
            Symbol.debug_print("_add_function_params:")
        # Note: we may be called from _fill_empty, so the symbols we want
        #       to add may actually already be present (as empty symbols).

        # add symbols for function parameters, if any
        if self.declaration is not None and self.declaration.function_params is not None:
            for p in self.declaration.function_params:
                if p.arg is None:
                    continue
                nn = p.arg.name
                if nn is None:
                    continue
                # (comparing to the template params: we have checked that we are a declaration)
                decl = ASTDeclaration('functionParam', None, p)
                assert not nn.rooted
                assert len(nn.names) == 1
                self._add_symbols(nn, decl, self.docname, self.line)
        if Symbol.debug_lookup:
            Symbol.debug_indent -= 1

    def remove(self) -> None:
        if self.parent:
            self.parent._remove_child(self)
            self.parent = None

    def clear_doc(self, docname: str) -> None:
        if docname not in self._children_by_docname:
            for child in self._children:
                child.clear_doc(docname)
            return

        children: dict[str, Symbol] = self._children_by_docname.pop(docname)
        for child in children.values():
            child.declaration = None
            child.docname = None
            child.line = None
            if child.siblingAbove is not None:
                child.siblingAbove.siblingBelow = child.siblingBelow
            if child.siblingBelow is not None:
                child.siblingBelow.siblingAbove = child.siblingAbove
            child.siblingAbove = None
            child.siblingBelow = None
            self._remove_child(child)

    def get_all_symbols(self) -> Iterator[Symbol]:
        yield self
        for sChild in self._children:
            yield from sChild.get_all_symbols()

    @property
    def children(self) -> Iterator[Symbol]:
        yield from self._children

    def get_lookup_key(self) -> LookupKey:
        # The pickle files for the environment and for each document are distinct.
        # The environment has all the symbols, but the documents has xrefs that
        # must know their scope. A lookup key is essentially a specification of
        # how to find a specific symbol.
        symbols = []
        s = self
        while s.parent:
            symbols.append(s)
            s = s.parent
        symbols.reverse()
        key = []
        for s in symbols:
            if s.declaration is not None:
                # TODO: do we need the ID?
                key.append((s.ident, s.declaration.get_newest_id()))
            else:
                key.append((s.ident, None))
        return LookupKey(key)

    def get_full_nested_name(self) -> ASTNestedName:
        symbols = []
        s = self
        while s.parent:
            symbols.append(s)
            s = s.parent
        symbols.reverse()
        names = [s.ident for s in symbols]
        return ASTNestedName(names, rooted=False)

    def _symbol_lookup(
        self,
        nestedName: ASTNestedName,
        onMissingQualifiedSymbol: Callable[[Symbol, ASTIdentifier], Symbol | None],
        ancestorLookupType: str | None,
        matchSelf: bool,
        recurseInAnon: bool,
        searchInSiblings: bool,
    ) -> SymbolLookupResult | None:
        # TODO: further simplification from C++ to C
        # ancestorLookupType: if not None, specifies the target type of the lookup
        if Symbol.debug_lookup:
            Symbol.debug_indent += 1
            Symbol.debug_print("_symbol_lookup:")
            Symbol.debug_indent += 1
            Symbol.debug_print("self:")
            logger.debug(self.to_string(Symbol.debug_indent + 1, addEndNewline=False))
            Symbol.debug_print("nestedName:        ", nestedName)
            Symbol.debug_print("ancestorLookupType:", ancestorLookupType)
            Symbol.debug_print("matchSelf:         ", matchSelf)
            Symbol.debug_print("recurseInAnon:     ", recurseInAnon)
            Symbol.debug_print("searchInSiblings:  ", searchInSiblings)

        names = nestedName.names

        # find the right starting point for lookup
        parentSymbol = self
        if nestedName.rooted:
            while parentSymbol.parent is not None:
                parentSymbol = parentSymbol.parent

        if ancestorLookupType is not None:
            # walk up until we find the first identifier
            firstName = names[0]
            while parentSymbol.parent:
                if firstName.name in parentSymbol._children_by_name:
                    break
                parentSymbol = parentSymbol.parent

        if Symbol.debug_lookup:
            Symbol.debug_print("starting point:")
            logger.debug(parentSymbol.to_string(Symbol.debug_indent + 1, addEndNewline=False))

        # and now the actual lookup
        for ident in names[:-1]:
            name = ident.name
            if name in parentSymbol._children_by_name:
                symbol = parentSymbol._children_by_name[name]
            else:
                symbol = onMissingQualifiedSymbol(parentSymbol, ident)
                if symbol is None:
                    if Symbol.debug_lookup:
                        Symbol.debug_indent -= 2
                    return None
            parentSymbol = symbol

        if Symbol.debug_lookup:
            Symbol.debug_print("handle last name from:")
            logger.debug(parentSymbol.to_string(Symbol.debug_indent + 1, addEndNewline=False))

        # handle the last name
        ident = names[-1]
        name = ident.name
        symbol = parentSymbol._children_by_name.get(name)
        if not symbol and recurseInAnon:
            for child in parentSymbol._anon_children:
                if name in child._children_by_name:
                    symbol = child._children_by_name[name]
                    break

        if Symbol.debug_lookup:
            Symbol.debug_indent -= 2

        result = [symbol] if symbol else []
        return SymbolLookupResult(result, parentSymbol, ident)

    def _add_symbols(
        self,
        nestedName: ASTNestedName,
        declaration: ASTDeclaration | None,
        docname: str | None,
        line: int | None,
    ) -> Symbol:
        # TODO: further simplification from C++ to C
        # Used for adding a whole path of symbols, where the last may or may not
        # be an actual declaration.

        if Symbol.debug_lookup:
            Symbol.debug_indent += 1
            Symbol.debug_print("_add_symbols:")
            Symbol.debug_indent += 1
            Symbol.debug_print("nn:       ", nestedName)
            Symbol.debug_print("decl:     ", declaration)
            Symbol.debug_print(f"location: {docname}:{line}")

        def onMissingQualifiedSymbol(parentSymbol: Symbol, ident: ASTIdentifier) -> Symbol:
            if Symbol.debug_lookup:
                Symbol.debug_indent += 1
                Symbol.debug_print("_add_symbols, onMissingQualifiedSymbol:")
                Symbol.debug_indent += 1
                Symbol.debug_print("ident: ", ident)
                Symbol.debug_indent -= 2
            return Symbol(parent=parentSymbol, ident=ident,
                          declaration=None, docname=None, line=None)

        lookupResult = self._symbol_lookup(nestedName,
                                           onMissingQualifiedSymbol,
                                           ancestorLookupType=None,
                                           matchSelf=False,
                                           recurseInAnon=False,
                                           searchInSiblings=False)
        assert lookupResult is not None  # we create symbols all the way, so that can't happen
        symbols = list(lookupResult.symbols)
        if len(symbols) == 0:
            if Symbol.debug_lookup:
                Symbol.debug_print("_add_symbols, result, no symbol:")
                Symbol.debug_indent += 1
                Symbol.debug_print("ident:       ", lookupResult.ident)
                Symbol.debug_print("declaration: ", declaration)
                Symbol.debug_print(f"location:    {docname}:{line}")
                Symbol.debug_indent -= 1
            symbol = Symbol(parent=lookupResult.parentSymbol,
                            ident=lookupResult.ident,
                            declaration=declaration,
                            docname=docname, line=line)
            if Symbol.debug_lookup:
                Symbol.debug_indent -= 2
            return symbol

        if Symbol.debug_lookup:
            Symbol.debug_print("_add_symbols, result, symbols:")
            Symbol.debug_indent += 1
            Symbol.debug_print("number symbols:", len(symbols))
            Symbol.debug_indent -= 1

        if not declaration:
            if Symbol.debug_lookup:
                Symbol.debug_print("no declaration")
                Symbol.debug_indent -= 2
            # good, just a scope creation
            # TODO: what if we have more than one symbol?
            return symbols[0]

        noDecl = []
        withDecl = []
        dupDecl = []
        for s in symbols:
            if s.declaration is None:
                noDecl.append(s)
            elif s.isRedeclaration:
                dupDecl.append(s)
            else:
                withDecl.append(s)
        if Symbol.debug_lookup:
            Symbol.debug_print("#noDecl:  ", len(noDecl))
            Symbol.debug_print("#withDecl:", len(withDecl))
            Symbol.debug_print("#dupDecl: ", len(dupDecl))

        # With partial builds we may start with a large symbol tree stripped of declarations.
        # Essentially any combination of noDecl, withDecl, and dupDecls seems possible.
        # TODO: make partial builds fully work. What should happen when the primary symbol gets
        #  deleted, and other duplicates exist? The full document should probably be rebuild.

        # First check if one of those with a declaration matches.
        # If it's a function, we need to compare IDs,
        # otherwise there should be only one symbol with a declaration.
        def makeCandSymbol() -> Symbol:
            if Symbol.debug_lookup:
                Symbol.debug_print("begin: creating candidate symbol")
            symbol = Symbol(parent=lookupResult.parentSymbol,
                            ident=lookupResult.ident,
                            declaration=declaration,
                            docname=docname, line=line)
            if Symbol.debug_lookup:
                Symbol.debug_print("end:   creating candidate symbol")
            return symbol

        if len(withDecl) == 0:
            candSymbol = None
        else:
            candSymbol = makeCandSymbol()

            def handleDuplicateDeclaration(symbol: Symbol, candSymbol: Symbol) -> None:
                if Symbol.debug_lookup:
                    Symbol.debug_indent += 1
                    Symbol.debug_print("redeclaration")
                    Symbol.debug_indent -= 1
                    Symbol.debug_indent -= 2
                # Redeclaration of the same symbol.
                # Let the new one be there, but raise an error to the client
                # so it can use the real symbol as subscope.
                # This will probably result in a duplicate id warning.
                candSymbol.isRedeclaration = True
                raise _DuplicateSymbolError(symbol, declaration)

            if declaration.objectType != "function":
                assert len(withDecl) <= 1
                handleDuplicateDeclaration(withDecl[0], candSymbol)
                # (not reachable)

            # a function, so compare IDs
            candId = declaration.get_newest_id()
            if Symbol.debug_lookup:
                Symbol.debug_print("candId:", candId)
            for symbol in withDecl:
                oldId = symbol.declaration.get_newest_id()
                if Symbol.debug_lookup:
                    Symbol.debug_print("oldId: ", oldId)
                if candId == oldId:
                    handleDuplicateDeclaration(symbol, candSymbol)
                    # (not reachable)
            # no candidate symbol found with matching ID
        # if there is an empty symbol, fill that one
        if len(noDecl) == 0:
            if Symbol.debug_lookup:
                Symbol.debug_print(
                    "no match, no empty, candSybmol is not None?:", candSymbol is not None,
                )
                Symbol.debug_indent -= 2
            if candSymbol is not None:
                return candSymbol
            else:
                return makeCandSymbol()
        else:
            if Symbol.debug_lookup:
                Symbol.debug_print(
                    "no match, but fill an empty declaration, candSybmol is not None?:",
                    candSymbol is not None)
                Symbol.debug_indent -= 2
            if candSymbol is not None:
                candSymbol.remove()
            # assert len(noDecl) == 1
            # TODO: enable assertion when we at some point find out how to do cleanup
            # for now, just take the first one, it should work fine ... right?
            symbol = noDecl[0]
            # If someone first opened the scope, and then later
            # declares it, e.g,
            # .. namespace:: Test
            # .. namespace:: nullptr
            # .. class:: Test
            symbol._fill_empty(declaration, docname, line)
            return symbol

    def merge_with(self, other: Symbol, docnames: list[str],
                   env: BuildEnvironment) -> None:
        if Symbol.debug_lookup:
            Symbol.debug_indent += 1
            Symbol.debug_print("merge_with:")

        assert other is not None
        for otherChild in other._children:
            otherName = otherChild.ident.name
            if otherName not in self._children_by_name:
                # TODO: hmm, should we prune by docnames?
                otherChild.parent = self
                self._add_child(otherChild)
                otherChild._assert_invariants()
                continue
            ourChild = self._children_by_name[otherName]
            if otherChild.declaration and otherChild.docname in docnames:
                if not ourChild.declaration:
                    ourChild._fill_empty(otherChild.declaration,
                                         otherChild.docname, otherChild.line)
                elif ourChild.docname != otherChild.docname:
                    name = str(ourChild.declaration)
                    msg = __("Duplicate C declaration, also defined at %s:%s.\n"
                             "Declaration is '.. c:%s:: %s'.")
                    msg = msg % (ourChild.docname, ourChild.line,
                                 ourChild.declaration.directiveType, name)
                    logger.warning(msg, location=(otherChild.docname, otherChild.line))
                else:
                    # Both have declarations, and in the same docname.
                    # This can apparently happen, it should be safe to
                    # just ignore it, right?
                    pass
            ourChild.merge_with(otherChild, docnames, env)

        if Symbol.debug_lookup:
            Symbol.debug_indent -= 1

    def add_name(self, nestedName: ASTNestedName) -> Symbol:
        if Symbol.debug_lookup:
            Symbol.debug_indent += 1
            Symbol.debug_print("add_name:")
        res = self._add_symbols(nestedName, declaration=None, docname=None, line=None)
        if Symbol.debug_lookup:
            Symbol.debug_indent -= 1
        return res

    def add_declaration(self, declaration: ASTDeclaration,
                        docname: str, line: int) -> Symbol:
        if Symbol.debug_lookup:
            Symbol.debug_indent += 1
            Symbol.debug_print("add_declaration:")
        assert declaration is not None
        assert docname is not None
        assert line is not None
        nestedName = declaration.name
        res = self._add_symbols(nestedName, declaration, docname, line)
        if Symbol.debug_lookup:
            Symbol.debug_indent -= 1
        return res

    def find_identifier(self, ident: ASTIdentifier,
                        matchSelf: bool, recurseInAnon: bool, searchInSiblings: bool,
                        ) -> Symbol | None:
        if Symbol.debug_lookup:
            Symbol.debug_indent += 1
            Symbol.debug_print("find_identifier:")
            Symbol.debug_indent += 1
            Symbol.debug_print("ident:           ", ident)
            Symbol.debug_print("matchSelf:       ", matchSelf)
            Symbol.debug_print("recurseInAnon:   ", recurseInAnon)
            Symbol.debug_print("searchInSiblings:", searchInSiblings)
            logger.debug(self.to_string(Symbol.debug_indent + 1, addEndNewline=False))
            Symbol.debug_indent -= 2
        current = self
        while current is not None:
            if Symbol.debug_lookup:
                Symbol.debug_indent += 2
                Symbol.debug_print("trying:")
                logger.debug(current.to_string(Symbol.debug_indent + 1, addEndNewline=False))
                Symbol.debug_indent -= 2
            if matchSelf and current.ident == ident:
                return current
            name = ident.name
            if name in current._children_by_name:
                return current._children_by_name[name]
            if recurseInAnon:
                for child in current._anon_children:
                    if name in child._children_by_name:
                        return child._children_by_name[name]
            if not searchInSiblings:
                break
            current = current.siblingAbove
        return None

    def direct_lookup(self, key: LookupKey) -> Symbol | None:
        if Symbol.debug_lookup:
            Symbol.debug_indent += 1
            Symbol.debug_print("direct_lookup:")
            Symbol.debug_indent += 1
        s = self
        for ident, id_ in key.data:
            s = s._children_by_name.get(ident.name)
            if Symbol.debug_lookup:
                Symbol.debug_print("name:          ", ident.name)
                Symbol.debug_print("id:            ", id_)
                if s is not None:
                    logger.debug(s.to_string(Symbol.debug_indent + 1, addEndNewline=False))
                else:
                    Symbol.debug_print("not found")
            if s is None:
                break
        if Symbol.debug_lookup:
            Symbol.debug_indent -= 2
        return s

    def find_declaration(self, nestedName: ASTNestedName, typ: str,
                         matchSelf: bool, recurseInAnon: bool) -> Symbol | None:
        # templateShorthand: missing template parameter lists for templates is ok
        if Symbol.debug_lookup:
            Symbol.debug_indent += 1
            Symbol.debug_print("find_declaration:")

        def onMissingQualifiedSymbol(
            parentSymbol: Symbol,
            ident: ASTIdentifier,
        ) -> Symbol | None:
            return None

        lookupResult = self._symbol_lookup(nestedName,
                                           onMissingQualifiedSymbol,
                                           ancestorLookupType=typ,
                                           matchSelf=matchSelf,
                                           recurseInAnon=recurseInAnon,
                                           searchInSiblings=False)
        if Symbol.debug_lookup:
            Symbol.debug_indent -= 1
        if lookupResult is None:
            return None

        symbols = list(lookupResult.symbols)
        if len(symbols) == 0:
            return None
        return symbols[0]

    def to_string(self, indent: int, *, addEndNewline: bool = True) -> str:
        res = [Symbol.debug_indent_string * indent]
        if not self.parent:
            res.append('::')
        else:
            if self.ident:
                res.append(self.ident.name)
            else:
                res.append(str(self.declaration))
            if self.declaration:
                res.append(": ")
                if self.isRedeclaration:
                    res.append('!!duplicate!! ')
                res.append(str(self.declaration))
        if self.docname:
            res.append('\t(')
            res.append(self.docname)
            res.append(')')
        if addEndNewline:
            res.append('\n')
        return ''.join(res)

    def dump(self, indent: int) -> str:
        return ''.join([self.to_string(indent), *(c.dump(indent + 1) for c in self._children)])
