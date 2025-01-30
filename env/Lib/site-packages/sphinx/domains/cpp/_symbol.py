from __future__ import annotations

from typing import TYPE_CHECKING, Any, NoReturn

from sphinx.domains.cpp._ast import (
    ASTDeclaration,
    ASTIdentifier,
    ASTNestedName,
    ASTNestedNameElement,
    ASTOperator,
    ASTTemplateArgs,
    ASTTemplateDeclarationPrefix,
    ASTTemplateIntroduction,
    ASTTemplateParams,
)
from sphinx.locale import __
from sphinx.util import logging

if TYPE_CHECKING:
    from collections.abc import Callable, Iterator

    from sphinx.environment import BuildEnvironment

logger = logging.getLogger(__name__)


class _DuplicateSymbolError(Exception):
    def __init__(self, symbol: Symbol, declaration: ASTDeclaration) -> None:
        assert symbol
        assert declaration
        self.symbol = symbol
        self.declaration = declaration

    def __str__(self) -> str:
        return "Internal C++ duplicate symbol error:\n%s" % self.symbol.dump(0)


class SymbolLookupResult:
    def __init__(self, symbols: Iterator[Symbol], parentSymbol: Symbol,
                 identOrOp: ASTIdentifier | ASTOperator, templateParams: Any,
                 templateArgs: ASTTemplateArgs) -> None:
        self.symbols = symbols
        self.parentSymbol = parentSymbol
        self.identOrOp = identOrOp
        self.templateParams = templateParams
        self.templateArgs = templateArgs


class LookupKey:
    def __init__(self, data: list[tuple[ASTNestedNameElement,
                                        ASTTemplateParams | ASTTemplateIntroduction,
                                        str]]) -> None:
        self.data = data


def _is_specialization(templateParams: ASTTemplateParams | ASTTemplateIntroduction,
                       templateArgs: ASTTemplateArgs) -> bool:
    # Checks if `templateArgs` does not exactly match `templateParams`.
    # the names of the template parameters must be given exactly as args
    # and params that are packs must in the args be the name expanded
    if len(templateParams.params) != len(templateArgs.args):
        return True
    # having no template params and no arguments is also a specialization
    if len(templateParams.params) == 0:
        return True
    for i in range(len(templateParams.params)):
        param = templateParams.params[i]
        arg = templateArgs.args[i]
        # TODO: doing this by string manipulation is probably not the most efficient
        paramName = str(param.name)
        argTxt = str(arg)
        isArgPackExpansion = argTxt.endswith('...')
        if param.isPack != isArgPackExpansion:
            return True
        argName = argTxt[:-3] if isArgPackExpansion else argTxt
        if paramName != argName:
            return True
    return False


class Symbol:
    debug_indent = 0
    debug_indent_string = "  "
    debug_lookup = False  # overridden by the corresponding config value
    debug_show_tree = False  # overridden by the corresponding config value

    def __copy__(self) -> NoReturn:
        raise AssertionError  # shouldn't happen

    def __deepcopy__(self, memo: Any) -> Symbol:
        if self.parent:
            raise AssertionError  # shouldn't happen
        # the domain base class makes a copy of the initial data, which is fine
        return Symbol(None, None, None, None, None, None, None)

    @staticmethod
    def debug_print(*args: Any) -> None:
        logger.debug(Symbol.debug_indent_string * Symbol.debug_indent, end="")
        logger.debug(*args)

    def _assert_invariants(self) -> None:
        if not self.parent:
            # parent == None means global scope, so declaration means a parent
            assert not self.identOrOp
            assert not self.templateParams
            assert not self.templateArgs
            assert not self.declaration
            assert not self.docname
        else:
            if self.declaration:
                assert self.docname

    def __setattr__(self, key: str, value: Any) -> None:
        if key == "children":
            raise AssertionError
        return super().__setattr__(key, value)

    def __init__(self, parent: Symbol | None,
                 identOrOp: ASTIdentifier | ASTOperator | None,
                 templateParams: ASTTemplateParams | ASTTemplateIntroduction | None,
                 templateArgs: Any, declaration: ASTDeclaration | None,
                 docname: str | None, line: int | None) -> None:
        self.parent = parent
        # declarations in a single directive are linked together
        self.siblingAbove: Symbol | None = None
        self.siblingBelow: Symbol | None = None
        self.identOrOp = identOrOp
        # Ensure the same symbol for `A` is created for:
        #
        #     .. cpp:class:: template <typename T> class A
        #
        # and
        #
        #     .. cpp:function:: template <typename T> int A<T>::foo()
        if (templateArgs is not None and
                not _is_specialization(templateParams, templateArgs)):
            templateArgs = None
        self.templateParams = templateParams  # template<templateParams>
        self.templateArgs = templateArgs  # identifier<templateArgs>
        self.declaration = declaration
        self.docname = docname
        self.line = line
        self.isRedeclaration = False
        self._assert_invariants()

        # Remember to modify Symbol.remove if modifications to the parent change.
        self._children: list[Symbol] = []
        self._anonChildren: list[Symbol] = []
        # note: _children includes _anonChildren
        if self.parent:
            self.parent._children.append(self)
        if self.declaration:
            self.declaration.symbol = self

        # Do symbol addition after self._children has been initialised.
        self._add_template_and_function_params()

    def __repr__(self) -> str:
        return f'<Symbol {self.to_string(indent=0)!r}>'

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
        self._add_template_and_function_params()

    def _add_template_and_function_params(self) -> None:
        if Symbol.debug_lookup:
            Symbol.debug_indent += 1
            Symbol.debug_print("_add_template_and_function_params:")
        # Note: we may be called from _fill_empty, so the symbols we want
        #       to add may actually already be present (as empty symbols).

        # add symbols for the template params
        if self.templateParams:
            for tp in self.templateParams.params:
                if not tp.get_identifier():
                    continue
                # only add a declaration if we our self are from a declaration
                if self.declaration:
                    decl = ASTDeclaration(objectType='templateParam', declaration=tp)
                else:
                    decl = None
                nne = ASTNestedNameElement(tp.get_identifier(), None)
                nn = ASTNestedName([nne], [False], rooted=False)
                self._add_symbols(nn, [], decl, self.docname, self.line)
        # add symbols for function parameters, if any
        if self.declaration is not None and self.declaration.function_params is not None:
            for fp in self.declaration.function_params:
                if fp.arg is None:
                    continue
                nn = fp.arg.name
                if nn is None:
                    continue
                # (comparing to the template params: we have checked that we are a declaration)
                decl = ASTDeclaration(objectType='functionParam', declaration=fp)
                assert not nn.rooted
                assert len(nn.names) == 1
                self._add_symbols(nn, [], decl, self.docname, self.line)
        if Symbol.debug_lookup:
            Symbol.debug_indent -= 1

    def remove(self) -> None:
        if self.parent is None:
            return
        assert self in self.parent._children
        self.parent._children.remove(self)
        self.parent = None

    def clear_doc(self, docname: str) -> None:
        newChildren: list[Symbol] = []
        for sChild in self._children:
            sChild.clear_doc(docname)
            if sChild.declaration and sChild.docname == docname:
                sChild.declaration = None
                sChild.docname = None
                sChild.line = None
                if sChild.siblingAbove is not None:
                    sChild.siblingAbove.siblingBelow = sChild.siblingBelow
                if sChild.siblingBelow is not None:
                    sChild.siblingBelow.siblingAbove = sChild.siblingAbove
                sChild.siblingAbove = None
                sChild.siblingBelow = None
            newChildren.append(sChild)
        self._children = newChildren

    def get_all_symbols(self) -> Iterator[Any]:
        yield self
        for sChild in self._children:
            yield from sChild.get_all_symbols()

    @property
    def children_recurse_anon(self) -> Iterator[Symbol]:
        for c in self._children:
            yield c
            if not c.identOrOp.is_anon():
                continue

            yield from c.children_recurse_anon

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
            nne = ASTNestedNameElement(s.identOrOp, s.templateArgs)
            if s.declaration is not None:
                key.append((nne, s.templateParams, s.declaration.get_newest_id()))
            else:
                key.append((nne, s.templateParams, None))
        return LookupKey(key)

    def get_full_nested_name(self) -> ASTNestedName:
        symbols = []
        s = self
        while s.parent:
            symbols.append(s)
            s = s.parent
        symbols.reverse()
        names = []
        templates = []
        for s in symbols:
            names.append(ASTNestedNameElement(s.identOrOp, s.templateArgs))
            templates.append(False)
        return ASTNestedName(names, templates, rooted=False)

    def _find_first_named_symbol(self, identOrOp: ASTIdentifier | ASTOperator,
                                 templateParams: ASTTemplateParams | ASTTemplateIntroduction,
                                 templateArgs: ASTTemplateArgs | None,
                                 templateShorthand: bool, matchSelf: bool,
                                 recurseInAnon: bool, correctPrimaryTemplateArgs: bool,
                                 ) -> Symbol | None:
        if Symbol.debug_lookup:
            Symbol.debug_print("_find_first_named_symbol ->")
        res = self._find_named_symbols(identOrOp, templateParams, templateArgs,
                                       templateShorthand, matchSelf, recurseInAnon,
                                       correctPrimaryTemplateArgs,
                                       searchInSiblings=False)
        try:
            return next(res)
        except StopIteration:
            return None

    def _find_named_symbols(self, identOrOp: ASTIdentifier | ASTOperator,
                            templateParams: ASTTemplateParams | ASTTemplateIntroduction,
                            templateArgs: ASTTemplateArgs,
                            templateShorthand: bool, matchSelf: bool,
                            recurseInAnon: bool, correctPrimaryTemplateArgs: bool,
                            searchInSiblings: bool) -> Iterator[Symbol]:
        if Symbol.debug_lookup:
            Symbol.debug_indent += 1
            Symbol.debug_print("_find_named_symbols:")
            Symbol.debug_indent += 1
            Symbol.debug_print("self:")
            logger.debug(self.to_string(Symbol.debug_indent + 1), end="")
            Symbol.debug_print("identOrOp:                  ", identOrOp)
            Symbol.debug_print("templateParams:             ", templateParams)
            Symbol.debug_print("templateArgs:               ", templateArgs)
            Symbol.debug_print("templateShorthand:          ", templateShorthand)
            Symbol.debug_print("matchSelf:                  ", matchSelf)
            Symbol.debug_print("recurseInAnon:              ", recurseInAnon)
            Symbol.debug_print("correctPrimaryTemplateAargs:", correctPrimaryTemplateArgs)
            Symbol.debug_print("searchInSiblings:           ", searchInSiblings)

        if correctPrimaryTemplateArgs:
            if templateParams is not None and templateArgs is not None:
                # If both are given, but it's not a specialization, then do lookup as if
                # there is no argument list.
                # For example: template<typename T> int A<T>::var;
                if not _is_specialization(templateParams, templateArgs):
                    templateArgs = None

        def matches(s: Symbol) -> bool:
            if s.identOrOp != identOrOp:
                return False
            if (s.templateParams is None) != (templateParams is None):
                if templateParams is not None:
                    # we query with params, they must match params
                    return False
                if not templateShorthand:
                    # we don't query with params, and we do care about them
                    return False
            if templateParams:
                # TODO: do better comparison
                if str(s.templateParams) != str(templateParams):
                    return False
            if (s.templateArgs is None) != (templateArgs is None):
                return False
            if s.templateArgs:
                # TODO: do better comparison
                if str(s.templateArgs) != str(templateArgs):
                    return False
            return True

        def candidates() -> Iterator[Symbol]:
            s = self
            if Symbol.debug_lookup:
                Symbol.debug_print("searching in self:")
                logger.debug(s.to_string(Symbol.debug_indent + 1), end="")
            while True:
                if matchSelf:
                    yield s
                if recurseInAnon:
                    yield from s.children_recurse_anon
                else:
                    yield from s._children

                if s.siblingAbove is None:
                    break
                s = s.siblingAbove
                if Symbol.debug_lookup:
                    Symbol.debug_print("searching in sibling:")
                    logger.debug(s.to_string(Symbol.debug_indent + 1), end="")

        for s in candidates():
            if Symbol.debug_lookup:
                Symbol.debug_print("candidate:")
                logger.debug(s.to_string(Symbol.debug_indent + 1), end="")
            if matches(s):
                if Symbol.debug_lookup:
                    Symbol.debug_indent += 1
                    Symbol.debug_print("matches")
                    Symbol.debug_indent -= 3
                yield s
                if Symbol.debug_lookup:
                    Symbol.debug_indent += 2
        if Symbol.debug_lookup:
            Symbol.debug_indent -= 2

    def _symbol_lookup(
        self,
        nestedName: ASTNestedName,
        templateDecls: list[Any],
        onMissingQualifiedSymbol: Callable[
            [Symbol, ASTIdentifier | ASTOperator, Any, ASTTemplateArgs], Symbol | None,
        ],
        strictTemplateParamArgLists: bool, ancestorLookupType: str,
        templateShorthand: bool, matchSelf: bool,
        recurseInAnon: bool, correctPrimaryTemplateArgs: bool,
        searchInSiblings: bool,
    ) -> SymbolLookupResult:
        # ancestorLookupType: if not None, specifies the target type of the lookup
        if Symbol.debug_lookup:
            Symbol.debug_indent += 1
            Symbol.debug_print("_symbol_lookup:")
            Symbol.debug_indent += 1
            Symbol.debug_print("self:")
            logger.debug(self.to_string(Symbol.debug_indent + 1), end="")
            Symbol.debug_print("nestedName:        ", nestedName)
            Symbol.debug_print("templateDecls:     ", ",".join(str(t) for t in templateDecls))
            Symbol.debug_print("strictTemplateParamArgLists:", strictTemplateParamArgLists)
            Symbol.debug_print("ancestorLookupType:", ancestorLookupType)
            Symbol.debug_print("templateShorthand: ", templateShorthand)
            Symbol.debug_print("matchSelf:         ", matchSelf)
            Symbol.debug_print("recurseInAnon:     ", recurseInAnon)
            Symbol.debug_print("correctPrimaryTemplateArgs: ", correctPrimaryTemplateArgs)
            Symbol.debug_print("searchInSiblings:  ", searchInSiblings)

        if strictTemplateParamArgLists:
            # Each template argument list must have a template parameter list.
            # But to declare a template there must be an additional template parameter list.
            assert (nestedName.num_templates() == len(templateDecls) or
                    nestedName.num_templates() + 1 == len(templateDecls))
        else:
            assert len(templateDecls) <= nestedName.num_templates() + 1

        names = nestedName.names

        # find the right starting point for lookup
        parentSymbol = self
        if nestedName.rooted:
            while parentSymbol.parent:
                parentSymbol = parentSymbol.parent
        if ancestorLookupType is not None:
            # walk up until we find the first identifier
            firstName = names[0]
            if not firstName.is_operator():
                while parentSymbol.parent:
                    if parentSymbol.find_identifier(firstName.identOrOp,
                                                    matchSelf=matchSelf,
                                                    recurseInAnon=recurseInAnon,
                                                    searchInSiblings=searchInSiblings):
                        # if we are in the scope of a constructor but wants to
                        # reference the class we need to walk one extra up
                        if (len(names) == 1 and ancestorLookupType == 'class' and matchSelf and
                                parentSymbol.parent and
                                parentSymbol.parent.identOrOp == firstName.identOrOp):
                            pass
                        else:
                            break
                    parentSymbol = parentSymbol.parent

        if Symbol.debug_lookup:
            Symbol.debug_print("starting point:")
            logger.debug(parentSymbol.to_string(Symbol.debug_indent + 1), end="")

        # and now the actual lookup
        iTemplateDecl = 0
        for name in names[:-1]:
            identOrOp = name.identOrOp
            templateArgs = name.templateArgs
            if strictTemplateParamArgLists:
                # there must be a parameter list
                if templateArgs:
                    assert iTemplateDecl < len(templateDecls)
                    templateParams = templateDecls[iTemplateDecl]
                    iTemplateDecl += 1
                else:
                    templateParams = None
            else:
                # take the next template parameter list if there is one
                # otherwise it's ok
                if templateArgs and iTemplateDecl < len(templateDecls):
                    templateParams = templateDecls[iTemplateDecl]
                    iTemplateDecl += 1
                else:
                    templateParams = None

            symbol = parentSymbol._find_first_named_symbol(
                identOrOp,
                templateParams, templateArgs,
                templateShorthand=templateShorthand,
                matchSelf=matchSelf,
                recurseInAnon=recurseInAnon,
                correctPrimaryTemplateArgs=correctPrimaryTemplateArgs)
            if symbol is None:
                symbol = onMissingQualifiedSymbol(parentSymbol, identOrOp,
                                                  templateParams, templateArgs)
                if symbol is None:
                    if Symbol.debug_lookup:
                        Symbol.debug_indent -= 2
                    return None
            # We have now matched part of a nested name, and need to match more
            # so even if we should matchSelf before, we definitely shouldn't
            # even more. (see also issue #2666)
            matchSelf = False
            parentSymbol = symbol

        if Symbol.debug_lookup:
            Symbol.debug_print("handle last name from:")
            logger.debug(parentSymbol.to_string(Symbol.debug_indent + 1), end="")

        # handle the last name
        name = names[-1]
        identOrOp = name.identOrOp
        templateArgs = name.templateArgs
        if iTemplateDecl < len(templateDecls):
            assert iTemplateDecl + 1 == len(templateDecls)
            templateParams = templateDecls[iTemplateDecl]
        else:
            assert iTemplateDecl == len(templateDecls)
            templateParams = None

        symbols = parentSymbol._find_named_symbols(
            identOrOp, templateParams, templateArgs,
            templateShorthand=templateShorthand, matchSelf=matchSelf,
            recurseInAnon=recurseInAnon, correctPrimaryTemplateArgs=False,
            searchInSiblings=searchInSiblings)
        if Symbol.debug_lookup:
            symbols = list(symbols)  # type: ignore[assignment]
            Symbol.debug_indent -= 2
        return SymbolLookupResult(symbols, parentSymbol,
                                  identOrOp, templateParams, templateArgs)

    def _add_symbols(
        self,
        nestedName: ASTNestedName,
        templateDecls: list[Any],
        declaration: ASTDeclaration | None,
        docname: str | None,
        line: int | None,
    ) -> Symbol:
        # Used for adding a whole path of symbols, where the last may or may not
        # be an actual declaration.

        if Symbol.debug_lookup:
            Symbol.debug_indent += 1
            Symbol.debug_print("_add_symbols:")
            Symbol.debug_indent += 1
            Symbol.debug_print("tdecls:", ",".join(str(t) for t in templateDecls))
            Symbol.debug_print("nn:       ", nestedName)
            Symbol.debug_print("decl:     ", declaration)
            Symbol.debug_print(f"location: {docname}:{line}")

        def onMissingQualifiedSymbol(parentSymbol: Symbol,
                                     identOrOp: ASTIdentifier | ASTOperator,
                                     templateParams: Any, templateArgs: ASTTemplateArgs,
                                     ) -> Symbol | None:
            if Symbol.debug_lookup:
                Symbol.debug_indent += 1
                Symbol.debug_print("_add_symbols, onMissingQualifiedSymbol:")
                Symbol.debug_indent += 1
                Symbol.debug_print("templateParams:", templateParams)
                Symbol.debug_print("identOrOp:     ", identOrOp)
                Symbol.debug_print("templateARgs:  ", templateArgs)
                Symbol.debug_indent -= 2
            return Symbol(parent=parentSymbol, identOrOp=identOrOp,
                          templateParams=templateParams,
                          templateArgs=templateArgs, declaration=None,
                          docname=None, line=None)

        lookupResult = self._symbol_lookup(nestedName, templateDecls,
                                           onMissingQualifiedSymbol,
                                           strictTemplateParamArgLists=True,
                                           ancestorLookupType=None,
                                           templateShorthand=False,
                                           matchSelf=False,
                                           recurseInAnon=False,
                                           correctPrimaryTemplateArgs=True,
                                           searchInSiblings=False)
        assert lookupResult is not None  # we create symbols all the way, so that can't happen
        symbols = list(lookupResult.symbols)
        if len(symbols) == 0:
            if Symbol.debug_lookup:
                Symbol.debug_print("_add_symbols, result, no symbol:")
                Symbol.debug_indent += 1
                Symbol.debug_print("templateParams:", lookupResult.templateParams)
                Symbol.debug_print("identOrOp:     ", lookupResult.identOrOp)
                Symbol.debug_print("templateArgs:  ", lookupResult.templateArgs)
                Symbol.debug_print("declaration:   ", declaration)
                Symbol.debug_print(f"location:      {docname}:{line}")
                Symbol.debug_indent -= 1
            symbol = Symbol(parent=lookupResult.parentSymbol,
                            identOrOp=lookupResult.identOrOp,
                            templateParams=lookupResult.templateParams,
                            templateArgs=lookupResult.templateArgs,
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
                            identOrOp=lookupResult.identOrOp,
                            templateParams=lookupResult.templateParams,
                            templateArgs=lookupResult.templateArgs,
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
                # but all existing must be functions as well,
                # otherwise we declare it to be a duplicate
                if symbol.declaration.objectType != 'function':
                    handleDuplicateDeclaration(symbol, candSymbol)
                    # (not reachable)
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
                Symbol.debug_print("no match, no empty")
                if candSymbol is not None:
                    Symbol.debug_print("result is already created candSymbol")
                else:
                    Symbol.debug_print("result is makeCandSymbol()")
                Symbol.debug_indent -= 2
            if candSymbol is not None:
                return candSymbol
            else:
                return makeCandSymbol()
        else:
            if Symbol.debug_lookup:
                Symbol.debug_print(
                    "no match, but fill an empty declaration, candSybmol is not None?:",
                    candSymbol is not None,
                )
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

        def unconditionalAdd(self: Symbol, otherChild: Symbol) -> None:
            # TODO: hmm, should we prune by docnames?
            self._children.append(otherChild)
            otherChild.parent = self
            otherChild._assert_invariants()

        if Symbol.debug_lookup:
            Symbol.debug_indent += 1
        for otherChild in other._children:
            if Symbol.debug_lookup:
                Symbol.debug_print("otherChild:\n", otherChild.to_string(Symbol.debug_indent))
                Symbol.debug_indent += 1
            if otherChild.isRedeclaration:
                unconditionalAdd(self, otherChild)
                if Symbol.debug_lookup:
                    Symbol.debug_print("isRedeclaration")
                    Symbol.debug_indent -= 1
                continue
            candiateIter = self._find_named_symbols(
                identOrOp=otherChild.identOrOp,
                templateParams=otherChild.templateParams,
                templateArgs=otherChild.templateArgs,
                templateShorthand=False, matchSelf=False,
                recurseInAnon=False, correctPrimaryTemplateArgs=False,
                searchInSiblings=False)
            candidates = list(candiateIter)

            if Symbol.debug_lookup:
                Symbol.debug_print("raw candidate symbols:", len(candidates))
            symbols = [s for s in candidates if not s.isRedeclaration]
            if Symbol.debug_lookup:
                Symbol.debug_print("non-duplicate candidate symbols:", len(symbols))

            if len(symbols) == 0:
                unconditionalAdd(self, otherChild)
                if Symbol.debug_lookup:
                    Symbol.debug_indent -= 1
                continue

            ourChild = None
            if otherChild.declaration is None:
                if Symbol.debug_lookup:
                    Symbol.debug_print("no declaration in other child")
                ourChild = symbols[0]
            else:
                queryId = otherChild.declaration.get_newest_id()
                if Symbol.debug_lookup:
                    Symbol.debug_print("queryId:  ", queryId)
                for symbol in symbols:
                    if symbol.declaration is None:
                        if Symbol.debug_lookup:
                            Symbol.debug_print("empty candidate")
                        # if in the end we have non-matching, but have an empty one,
                        # then just continue with that
                        ourChild = symbol
                        continue
                    candId = symbol.declaration.get_newest_id()
                    if Symbol.debug_lookup:
                        Symbol.debug_print("candidate:", candId)
                    if candId == queryId:
                        ourChild = symbol
                        break
            if Symbol.debug_lookup:
                Symbol.debug_indent -= 1
            if ourChild is None:
                unconditionalAdd(self, otherChild)
                continue
            if otherChild.declaration and otherChild.docname in docnames:
                if not ourChild.declaration:
                    ourChild._fill_empty(otherChild.declaration,
                                         otherChild.docname, otherChild.line)
                elif ourChild.docname != otherChild.docname:
                    name = str(ourChild.declaration)
                    msg = __("Duplicate C++ declaration, also defined at %s:%s.\n"
                             "Declaration is '.. cpp:%s:: %s'.")
                    msg = msg % (ourChild.docname, ourChild.line,
                                 ourChild.declaration.directiveType, name)
                    logger.warning(msg, location=(otherChild.docname, otherChild.line))
                else:
                    if (otherChild.declaration.objectType ==
                            ourChild.declaration.objectType and
                            otherChild.declaration.objectType in
                            ('templateParam', 'functionParam') and
                            ourChild.parent.declaration == otherChild.parent.declaration):
                        # `ourChild` was just created during merging by the call
                        # to `_fill_empty` on the parent and can be ignored.
                        pass
                    else:
                        # Both have declarations, and in the same docname.
                        # This can apparently happen, it should be safe to
                        # just ignore it, right?
                        # Hmm, only on duplicate declarations, right?
                        msg = "Internal C++ domain error during symbol merging.\n"
                        msg += "ourChild:\n" + ourChild.to_string(1)
                        msg += "\notherChild:\n" + otherChild.to_string(1)
                        logger.warning(msg, location=otherChild.docname)
            ourChild.merge_with(otherChild, docnames, env)
        if Symbol.debug_lookup:
            Symbol.debug_indent -= 2

    def add_name(self, nestedName: ASTNestedName,
                 templatePrefix: ASTTemplateDeclarationPrefix | None = None) -> Symbol:
        if Symbol.debug_lookup:
            Symbol.debug_indent += 1
            Symbol.debug_print("add_name:")
        if templatePrefix:
            templateDecls = templatePrefix.templates
        else:
            templateDecls = []
        res = self._add_symbols(nestedName, templateDecls,
                                declaration=None, docname=None, line=None)
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
        if declaration.templatePrefix:
            templateDecls = declaration.templatePrefix.templates
        else:
            templateDecls = []
        res = self._add_symbols(nestedName, templateDecls, declaration, docname, line)
        if Symbol.debug_lookup:
            Symbol.debug_indent -= 1
        return res

    def find_identifier(self, identOrOp: ASTIdentifier | ASTOperator,
                        matchSelf: bool, recurseInAnon: bool, searchInSiblings: bool,
                        ) -> Symbol | None:
        if Symbol.debug_lookup:
            Symbol.debug_indent += 1
            Symbol.debug_print("find_identifier:")
            Symbol.debug_indent += 1
            Symbol.debug_print("identOrOp:       ", identOrOp)
            Symbol.debug_print("matchSelf:       ", matchSelf)
            Symbol.debug_print("recurseInAnon:   ", recurseInAnon)
            Symbol.debug_print("searchInSiblings:", searchInSiblings)
            logger.debug(self.to_string(Symbol.debug_indent + 1), end="")
            Symbol.debug_indent -= 2
        current = self
        while current is not None:
            if Symbol.debug_lookup:
                Symbol.debug_indent += 2
                Symbol.debug_print("trying:")
                logger.debug(current.to_string(Symbol.debug_indent + 1), end="")
                Symbol.debug_indent -= 2
            if matchSelf and current.identOrOp == identOrOp:
                return current
            children = current.children_recurse_anon if recurseInAnon else current._children
            for s in children:
                if s.identOrOp == identOrOp:
                    return s
            if not searchInSiblings:
                break
            current = current.siblingAbove
        return None

    def direct_lookup(self, key: LookupKey) -> Symbol:
        if Symbol.debug_lookup:
            Symbol.debug_indent += 1
            Symbol.debug_print("direct_lookup:")
            Symbol.debug_indent += 1
        s = self
        for name, templateParams, id_ in key.data:
            if id_ is not None:
                res = None
                for cand in s._children:
                    if cand.declaration is None:
                        continue
                    if cand.declaration.get_newest_id() == id_:
                        res = cand
                        break
                s = res
            else:
                identOrOp = name.identOrOp
                templateArgs = name.templateArgs
                s = s._find_first_named_symbol(identOrOp,
                                               templateParams, templateArgs,
                                               templateShorthand=False,
                                               matchSelf=False,
                                               recurseInAnon=False,
                                               correctPrimaryTemplateArgs=False)
            if Symbol.debug_lookup:
                Symbol.debug_print("name:          ", name)
                Symbol.debug_print("templateParams:", templateParams)
                Symbol.debug_print("id:            ", id_)
                if s is not None:
                    logger.debug(s.to_string(Symbol.debug_indent + 1), end="")
                else:
                    Symbol.debug_print("not found")
            if s is None:
                if Symbol.debug_lookup:
                    Symbol.debug_indent -= 2
                return None
        if Symbol.debug_lookup:
            Symbol.debug_indent -= 2
        return s

    def find_name(
        self,
        nestedName: ASTNestedName,
        templateDecls: list[Any],
        typ: str,
        templateShorthand: bool,
        matchSelf: bool,
        recurseInAnon: bool,
        searchInSiblings: bool,
    ) -> tuple[list[Symbol] | None, str]:
        # templateShorthand: missing template parameter lists for templates is ok
        # If the first component is None,
        # then the second component _may_ be a string explaining why.
        if Symbol.debug_lookup:
            Symbol.debug_indent += 1
            Symbol.debug_print("find_name:")
            Symbol.debug_indent += 1
            Symbol.debug_print("self:")
            logger.debug(self.to_string(Symbol.debug_indent + 1), end="")
            Symbol.debug_print("nestedName:       ", nestedName)
            Symbol.debug_print("templateDecls:    ", templateDecls)
            Symbol.debug_print("typ:              ", typ)
            Symbol.debug_print("templateShorthand:", templateShorthand)
            Symbol.debug_print("matchSelf:        ", matchSelf)
            Symbol.debug_print("recurseInAnon:    ", recurseInAnon)
            Symbol.debug_print("searchInSiblings: ", searchInSiblings)

        class QualifiedSymbolIsTemplateParam(Exception):
            pass

        def onMissingQualifiedSymbol(parentSymbol: Symbol,
                                     identOrOp: ASTIdentifier | ASTOperator,
                                     templateParams: Any,
                                     templateArgs: ASTTemplateArgs) -> Symbol | None:
            # TODO: Maybe search without template args?
            #       Though, the correctPrimaryTemplateArgs does
            #       that for primary templates.
            #       Is there another case where it would be good?
            if parentSymbol.declaration is not None:
                if parentSymbol.declaration.objectType == 'templateParam':
                    raise QualifiedSymbolIsTemplateParam
            return None

        try:
            lookupResult = self._symbol_lookup(nestedName, templateDecls,
                                               onMissingQualifiedSymbol,
                                               strictTemplateParamArgLists=False,
                                               ancestorLookupType=typ,
                                               templateShorthand=templateShorthand,
                                               matchSelf=matchSelf,
                                               recurseInAnon=recurseInAnon,
                                               correctPrimaryTemplateArgs=False,
                                               searchInSiblings=searchInSiblings)
        except QualifiedSymbolIsTemplateParam:
            return None, "templateParamInQualified"

        if lookupResult is None:
            # if it was a part of the qualification that could not be found
            if Symbol.debug_lookup:
                Symbol.debug_indent -= 2
            return None, None

        res = list(lookupResult.symbols)
        if len(res) != 0:
            if Symbol.debug_lookup:
                Symbol.debug_indent -= 2
            return res, None

        if lookupResult.parentSymbol.declaration is not None:
            if lookupResult.parentSymbol.declaration.objectType == 'templateParam':
                return None, "templateParamInQualified"

        # try without template params and args
        symbol = lookupResult.parentSymbol._find_first_named_symbol(
            lookupResult.identOrOp, None, None,
            templateShorthand=templateShorthand, matchSelf=matchSelf,
            recurseInAnon=recurseInAnon, correctPrimaryTemplateArgs=False)
        if Symbol.debug_lookup:
            Symbol.debug_indent -= 2
        if symbol is not None:
            return [symbol], None
        else:
            return None, None

    def find_declaration(self, declaration: ASTDeclaration, typ: str, templateShorthand: bool,
                         matchSelf: bool, recurseInAnon: bool) -> Symbol | None:
        # templateShorthand: missing template parameter lists for templates is ok
        if Symbol.debug_lookup:
            Symbol.debug_indent += 1
            Symbol.debug_print("find_declaration:")
        nestedName = declaration.name
        if declaration.templatePrefix:
            templateDecls = declaration.templatePrefix.templates
        else:
            templateDecls = []

        def onMissingQualifiedSymbol(parentSymbol: Symbol,
                                     identOrOp: ASTIdentifier | ASTOperator,
                                     templateParams: Any,
                                     templateArgs: ASTTemplateArgs) -> Symbol | None:
            return None

        lookupResult = self._symbol_lookup(nestedName, templateDecls,
                                           onMissingQualifiedSymbol,
                                           strictTemplateParamArgLists=False,
                                           ancestorLookupType=typ,
                                           templateShorthand=templateShorthand,
                                           matchSelf=matchSelf,
                                           recurseInAnon=recurseInAnon,
                                           correctPrimaryTemplateArgs=False,
                                           searchInSiblings=False)
        if Symbol.debug_lookup:
            Symbol.debug_indent -= 1
        if lookupResult is None:
            return None

        symbols = list(lookupResult.symbols)
        if len(symbols) == 0:
            return None

        querySymbol = Symbol(parent=lookupResult.parentSymbol,
                             identOrOp=lookupResult.identOrOp,
                             templateParams=lookupResult.templateParams,
                             templateArgs=lookupResult.templateArgs,
                             declaration=declaration,
                             docname='fakeDocnameForQuery',
                             line=42)
        queryId = declaration.get_newest_id()
        for symbol in symbols:
            if symbol.declaration is None:
                continue
            candId = symbol.declaration.get_newest_id()
            if candId == queryId:
                querySymbol.remove()
                return symbol
        querySymbol.remove()
        return None

    def to_string(self, indent: int) -> str:
        res = [Symbol.debug_indent_string * indent]
        if not self.parent:
            res.append('::')
        else:
            if self.templateParams:
                res.append(str(self.templateParams))
                res.append('\n')
                res.append(Symbol.debug_indent_string * indent)
            if self.identOrOp:
                res.append(str(self.identOrOp))
            else:
                res.append(str(self.declaration))
            if self.templateArgs:
                res.append(str(self.templateArgs))
            if self.declaration:
                res.append(": ")
                if self.isRedeclaration:
                    res.append('!!duplicate!! ')
                res.append("{" + self.declaration.objectType + "} ")
                res.append(str(self.declaration))
        if self.docname:
            res.append('\t(')
            res.append(self.docname)
            res.append(')')
        res.append('\n')
        return ''.join(res)

    def dump(self, indent: int) -> str:
        return ''.join([
            self.to_string(indent),
            *(c.dump(indent + 1) for c in self._children),
        ])
