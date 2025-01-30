from __future__ import annotations

import re
from typing import TYPE_CHECKING, Any

from sphinx.domains.cpp._ast import (
    ASTAlignofExpr,
    ASTArray,
    ASTAssignmentExpr,
    ASTBaseClass,
    ASTBinOpExpr,
    ASTBooleanLiteral,
    ASTBracedInitList,
    ASTCastExpr,
    ASTCharLiteral,
    ASTClass,
    ASTCommaExpr,
    ASTConcept,
    ASTConditionalExpr,
    ASTDeclaration,
    ASTDeclarator,
    ASTDeclaratorMemPtr,
    ASTDeclaratorNameBitField,
    ASTDeclaratorNameParamQual,
    ASTDeclaratorParamPack,
    ASTDeclaratorParen,
    ASTDeclaratorPtr,
    ASTDeclaratorRef,
    ASTDeclSpecs,
    ASTDeclSpecsSimple,
    ASTDeleteExpr,
    ASTEnum,
    ASTEnumerator,
    ASTExplicitCast,
    ASTExplicitSpec,
    ASTExpression,
    ASTFallbackExpr,
    ASTFoldExpr,
    ASTFunctionParameter,
    ASTIdentifier,
    ASTIdExpression,
    ASTInitializer,
    ASTLiteral,
    ASTNamespace,
    ASTNestedName,
    ASTNestedNameElement,
    ASTNewExpr,
    ASTNoexceptExpr,
    ASTNoexceptSpec,
    ASTNumberLiteral,
    ASTOperator,
    ASTOperatorBuildIn,
    ASTOperatorLiteral,
    ASTOperatorType,
    ASTPackExpansionExpr,
    ASTParametersQualifiers,
    ASTParenExpr,
    ASTParenExprList,
    ASTPointerLiteral,
    ASTPostfixArray,
    ASTPostfixCallExpr,
    ASTPostfixDec,
    ASTPostfixExpr,
    ASTPostfixInc,
    ASTPostfixMember,
    ASTPostfixMemberOfPointer,
    ASTPostfixOp,
    ASTRequiresClause,
    ASTSizeofExpr,
    ASTSizeofParamPack,
    ASTSizeofType,
    ASTStringLiteral,
    ASTTemplateArgConstant,
    ASTTemplateArgs,
    ASTTemplateDeclarationPrefix,
    ASTTemplateIntroduction,
    ASTTemplateIntroductionParameter,
    ASTTemplateKeyParamPackIdDefault,
    ASTTemplateParam,
    ASTTemplateParamConstrainedTypeWithInit,
    ASTTemplateParamNonType,
    ASTTemplateParams,
    ASTTemplateParamTemplateType,
    ASTTemplateParamType,
    ASTThisLiteral,
    ASTTrailingTypeSpec,
    ASTTrailingTypeSpecDecltype,
    ASTTrailingTypeSpecDecltypeAuto,
    ASTTrailingTypeSpecFundamental,
    ASTTrailingTypeSpecName,
    ASTType,
    ASTTypeId,
    ASTTypeUsing,
    ASTTypeWithInit,
    ASTUnaryOpExpr,
    ASTUnion,
    ASTUserDefinedLiteral,
)
from sphinx.domains.cpp._ids import (
    _expression_assignment_ops,
    _expression_bin_ops,
    _expression_unary_ops,
    _fold_operator_re,
    _id_explicit_cast,
    _keywords,
    _operator_re,
    _simple_type_specifiers_re,
    _string_re,
    _visibility_re,
    udl_identifier_re,
)
from sphinx.util import logging
from sphinx.util.cfamily import (
    ASTAttributeList,
    BaseParser,
    DefinitionError,
    UnsupportedMultiCharacterCharLiteral,
    binary_literal_re,
    char_literal_re,
    float_literal_re,
    float_literal_suffix_re,
    hex_literal_re,
    identifier_re,
    integer_literal_re,
    integers_literal_suffix_re,
    octal_literal_re,
)

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

logger = logging.getLogger(__name__)


class DefinitionParser(BaseParser):
    @property
    def language(self) -> str:
        return 'C++'

    @property
    def id_attributes(self) -> Sequence[str]:
        return self.config.cpp_id_attributes

    @property
    def paren_attributes(self) -> Sequence[str]:
        return self.config.cpp_paren_attributes

    def _parse_string(self) -> str:
        if self.current_char != '"':
            return None
        startPos = self.pos
        self.pos += 1
        escape = False
        while True:
            if self.eof:
                self.fail("Unexpected end during inside string.")
            elif self.current_char == '"' and not escape:
                self.pos += 1
                break
            elif self.current_char == '\\':
                escape = True
            else:
                escape = False
            self.pos += 1
        return self.definition[startPos:self.pos]

    def _parse_literal(self) -> ASTLiteral:
        # -> integer-literal
        #  | character-literal
        #  | floating-literal
        #  | string-literal
        #  | boolean-literal -> "false" | "true"
        #  | pointer-literal -> "nullptr"
        #  | user-defined-literal

        def _udl(literal: ASTLiteral) -> ASTLiteral:
            if not self.match(udl_identifier_re):
                return literal
            # hmm, should we care if it's a keyword?
            # it looks like GCC does not disallow keywords
            ident = ASTIdentifier(self.matched_text)
            return ASTUserDefinedLiteral(literal, ident)

        self.skip_ws()
        if self.skip_word('nullptr'):
            return ASTPointerLiteral()
        if self.skip_word('true'):
            return ASTBooleanLiteral(True)
        if self.skip_word('false'):
            return ASTBooleanLiteral(False)
        pos = self.pos
        if self.match(float_literal_re):
            hasSuffix = self.match(float_literal_suffix_re)
            floatLit = ASTNumberLiteral(self.definition[pos:self.pos])
            if hasSuffix:
                return floatLit
            else:
                return _udl(floatLit)
        for regex in (binary_literal_re, hex_literal_re,
                      integer_literal_re, octal_literal_re):
            if self.match(regex):
                hasSuffix = self.match(integers_literal_suffix_re)
                intLit = ASTNumberLiteral(self.definition[pos:self.pos])
                if hasSuffix:
                    return intLit
                else:
                    return _udl(intLit)

        string = self._parse_string()
        if string is not None:
            return _udl(ASTStringLiteral(string))

        # character-literal
        if self.match(char_literal_re):
            prefix = self.last_match.group(1)  # may be None when no prefix
            data = self.last_match.group(2)
            try:
                charLit = ASTCharLiteral(prefix, data)
            except UnicodeDecodeError as e:
                self.fail("Can not handle character literal. Internal error was: %s" % e)
            except UnsupportedMultiCharacterCharLiteral:
                self.fail("Can not handle character literal"
                          " resulting in multiple decoded characters.")
            return _udl(charLit)
        return None

    def _parse_fold_or_paren_expression(self) -> ASTExpression | None:
        # "(" expression ")"
        # fold-expression
        # -> ( cast-expression fold-operator ... )
        #  | ( ... fold-operator cast-expression )
        #  | ( cast-expression fold-operator ... fold-operator cast-expression
        if self.current_char != '(':
            return None
        self.pos += 1
        self.skip_ws()
        if self.skip_string_and_ws("..."):
            # ( ... fold-operator cast-expression )
            if not self.match(_fold_operator_re):
                self.fail("Expected fold operator after '...' in fold expression.")
            op = self.matched_text
            rightExpr = self._parse_cast_expression()
            if not self.skip_string(')'):
                self.fail("Expected ')' in end of fold expression.")
            return ASTFoldExpr(None, op, rightExpr)
        # try first parsing a unary right fold, or a binary fold
        pos = self.pos
        try:
            self.skip_ws()
            leftExpr = self._parse_cast_expression()
            self.skip_ws()
            if not self.match(_fold_operator_re):
                self.fail("Expected fold operator after left expression in fold expression.")
            op = self.matched_text
            self.skip_ws()
            if not self.skip_string_and_ws('...'):
                self.fail("Expected '...' after fold operator in fold expression.")
        except DefinitionError as eFold:
            self.pos = pos
            # fall back to a paren expression
            try:
                res = self._parse_expression()
                self.skip_ws()
                if not self.skip_string(')'):
                    self.fail("Expected ')' in end of parenthesized expression.")
            except DefinitionError as eExpr:
                raise self._make_multi_error([
                    (eFold, "If fold expression"),
                    (eExpr, "If parenthesized expression"),
                ], "Error in fold expression or parenthesized expression.") from eExpr
            return ASTParenExpr(res)
        # now it definitely is a fold expression
        if self.skip_string(')'):
            return ASTFoldExpr(leftExpr, op, None)
        if not self.match(_fold_operator_re):
            self.fail("Expected fold operator or ')' after '...' in fold expression.")
        if op != self.matched_text:
            self.fail("Operators are different in binary fold: '%s' and '%s'."
                      % (op, self.matched_text))
        rightExpr = self._parse_cast_expression()
        self.skip_ws()
        if not self.skip_string(')'):
            self.fail("Expected ')' to end binary fold expression.")
        return ASTFoldExpr(leftExpr, op, rightExpr)

    def _parse_primary_expression(self) -> ASTExpression:
        # literal
        # "this"
        # lambda-expression
        # "(" expression ")"
        # fold-expression
        # id-expression -> we parse this with _parse_nested_name
        self.skip_ws()
        res: ASTExpression = self._parse_literal()
        if res is not None:
            return res
        self.skip_ws()
        if self.skip_word("this"):
            return ASTThisLiteral()
        # TODO: try lambda expression
        res = self._parse_fold_or_paren_expression()
        if res is not None:
            return res
        nn = self._parse_nested_name()
        if nn is not None:
            return ASTIdExpression(nn)
        return None

    def _parse_initializer_list(self, name: str, open: str, close: str,
                                ) -> tuple[list[ASTExpression | ASTBracedInitList],
                                           bool]:
        # Parse open and close with the actual initializer-list in between
        # -> initializer-clause '...'[opt]
        #  | initializer-list ',' initializer-clause '...'[opt]
        self.skip_ws()
        if not self.skip_string_and_ws(open):
            return None, None
        if self.skip_string(close):
            return [], False

        exprs: list[ASTExpression | ASTBracedInitList] = []
        trailingComma = False
        while True:
            self.skip_ws()
            expr = self._parse_initializer_clause()
            self.skip_ws()
            if self.skip_string('...'):
                exprs.append(ASTPackExpansionExpr(expr))
            else:
                exprs.append(expr)
            self.skip_ws()
            if self.skip_string(close):
                break
            if not self.skip_string_and_ws(','):
                self.fail(f"Error in {name}, expected ',' or '{close}'.")
            if self.current_char == close == '}':
                self.pos += 1
                trailingComma = True
                break
        return exprs, trailingComma

    def _parse_paren_expression_list(self) -> ASTParenExprList:
        # -> '(' expression-list ')'
        # though, we relax it to also allow empty parens
        # as it's needed in some cases
        #
        # expression-list
        # -> initializer-list
        exprs, trailingComma = self._parse_initializer_list("parenthesized expression-list",
                                                            '(', ')')
        if exprs is None:
            return None
        return ASTParenExprList(exprs)

    def _parse_initializer_clause(self) -> ASTExpression | ASTBracedInitList:
        bracedInitList = self._parse_braced_init_list()
        if bracedInitList is not None:
            return bracedInitList
        return self._parse_assignment_expression(inTemplate=False)

    def _parse_braced_init_list(self) -> ASTBracedInitList:
        # -> '{' initializer-list ','[opt] '}'
        #  | '{' '}'
        exprs, trailingComma = self._parse_initializer_list("braced-init-list", '{', '}')
        if exprs is None:
            return None
        return ASTBracedInitList(exprs, trailingComma)

    def _parse_expression_list_or_braced_init_list(
        self,
    ) -> ASTParenExprList | ASTBracedInitList:
        paren = self._parse_paren_expression_list()
        if paren is not None:
            return paren
        return self._parse_braced_init_list()

    def _parse_postfix_expression(self) -> ASTPostfixExpr:
        # -> primary
        #  | postfix "[" expression "]"
        #  | postfix "[" braced-init-list [opt] "]"
        #  | postfix "(" expression-list [opt] ")"
        #  | postfix "." "template" [opt] id-expression
        #  | postfix "->" "template" [opt] id-expression
        #  | postfix "." pseudo-destructor-name
        #  | postfix "->" pseudo-destructor-name
        #  | postfix "++"
        #  | postfix "--"
        #  | simple-type-specifier "(" expression-list [opt] ")"
        #  | simple-type-specifier braced-init-list
        #  | typename-specifier "(" expression-list [opt] ")"
        #  | typename-specifier braced-init-list
        #  | "dynamic_cast" "<" type-id ">" "(" expression ")"
        #  | "static_cast" "<" type-id ">" "(" expression ")"
        #  | "reinterpret_cast" "<" type-id ">" "(" expression ")"
        #  | "const_cast" "<" type-id ">" "(" expression ")"
        #  | "typeid" "(" expression ")"
        #  | "typeid" "(" type-id ")"

        prefixType = None
        prefix: Any = None
        self.skip_ws()

        cast = None
        for c in _id_explicit_cast:
            if self.skip_word_and_ws(c):
                cast = c
                break
        if cast is not None:
            prefixType = "cast"
            if not self.skip_string("<"):
                self.fail("Expected '<' after '%s'." % cast)
            typ = self._parse_type(False)
            self.skip_ws()
            if not self.skip_string_and_ws(">"):
                self.fail("Expected '>' after type in '%s'." % cast)
            if not self.skip_string("("):
                self.fail("Expected '(' in '%s'." % cast)

            def parser() -> ASTExpression:
                return self._parse_expression()
            expr = self._parse_expression_fallback([')'], parser)
            self.skip_ws()
            if not self.skip_string(")"):
                self.fail("Expected ')' to end '%s'." % cast)
            prefix = ASTExplicitCast(cast, typ, expr)
        elif self.skip_word_and_ws("typeid"):
            prefixType = "typeid"
            if not self.skip_string_and_ws('('):
                self.fail("Expected '(' after 'typeid'.")
            pos = self.pos
            try:
                typ = self._parse_type(False)
                prefix = ASTTypeId(typ, isType=True)
                if not self.skip_string(')'):
                    self.fail("Expected ')' to end 'typeid' of type.")
            except DefinitionError as eType:
                self.pos = pos
                try:

                    def parser() -> ASTExpression:
                        return self._parse_expression()
                    expr = self._parse_expression_fallback([')'], parser)
                    prefix = ASTTypeId(expr, isType=False)
                    if not self.skip_string(')'):
                        self.fail("Expected ')' to end 'typeid' of expression.")
                except DefinitionError as eExpr:
                    self.pos = pos
                    header = "Error in 'typeid(...)'."
                    header += " Expected type or expression."
                    errors = []
                    errors.append((eType, "If type"))
                    errors.append((eExpr, "If expression"))
                    raise self._make_multi_error(errors, header) from eExpr
        else:  # a primary expression or a type
            pos = self.pos
            try:
                prefix = self._parse_primary_expression()
                prefixType = 'expr'
            except DefinitionError as eOuter:
                self.pos = pos
                try:
                    # we are potentially casting, so save parens for us
                    # TODO: hmm, would we need to try both with operatorCast and with None?
                    prefix = self._parse_type(False, 'operatorCast')
                    prefixType = 'typeOperatorCast'
                    #  | simple-type-specifier "(" expression-list [opt] ")"
                    #  | simple-type-specifier braced-init-list
                    #  | typename-specifier "(" expression-list [opt] ")"
                    #  | typename-specifier braced-init-list
                    self.skip_ws()
                    if self.current_char != '(' and self.current_char != '{':
                        self.fail("Expecting '(' or '{' after type in cast expression.")
                except DefinitionError as eInner:
                    self.pos = pos
                    header = "Error in postfix expression,"
                    header += " expected primary expression or type."
                    errors = []
                    errors.append((eOuter, "If primary expression"))
                    errors.append((eInner, "If type"))
                    raise self._make_multi_error(errors, header) from eInner

        # and now parse postfixes
        postFixes: list[ASTPostfixOp] = []
        while True:
            self.skip_ws()
            if prefixType in ('expr', 'cast', 'typeid'):
                if self.skip_string_and_ws('['):
                    expr = self._parse_expression()
                    self.skip_ws()
                    if not self.skip_string(']'):
                        self.fail("Expected ']' in end of postfix expression.")
                    postFixes.append(ASTPostfixArray(expr))
                    continue
                if self.skip_string('.'):
                    if self.skip_string('*'):
                        # don't steal the dot
                        self.pos -= 2
                    elif self.skip_string('..'):
                        # don't steal the dot
                        self.pos -= 3
                    else:
                        name = self._parse_nested_name()
                        postFixes.append(ASTPostfixMember(name))
                        continue
                if self.skip_string('->'):
                    if self.skip_string('*'):
                        # don't steal the arrow
                        self.pos -= 3
                    else:
                        name = self._parse_nested_name()
                        postFixes.append(ASTPostfixMemberOfPointer(name))
                        continue
                if self.skip_string('++'):
                    postFixes.append(ASTPostfixInc())
                    continue
                if self.skip_string('--'):
                    postFixes.append(ASTPostfixDec())
                    continue
            lst = self._parse_expression_list_or_braced_init_list()
            if lst is not None:
                postFixes.append(ASTPostfixCallExpr(lst))
                continue
            break
        return ASTPostfixExpr(prefix, postFixes)

    def _parse_unary_expression(self) -> ASTExpression:
        # -> postfix
        #  | "++" cast
        #  | "--" cast
        #  | unary-operator cast -> (* | & | + | - | ! | ~) cast
        # The rest:
        #  | "sizeof" unary
        #  | "sizeof" "(" type-id ")"
        #  | "sizeof" "..." "(" identifier ")"
        #  | "alignof" "(" type-id ")"
        #  | noexcept-expression -> noexcept "(" expression ")"
        #  | new-expression
        #  | delete-expression
        self.skip_ws()
        for op in _expression_unary_ops:
            # TODO: hmm, should we be able to backtrack here?
            if op[0] in 'cn':
                res = self.skip_word(op)
            else:
                res = self.skip_string(op)
            if res:
                expr = self._parse_cast_expression()
                return ASTUnaryOpExpr(op, expr)
        if self.skip_word_and_ws('sizeof'):
            if self.skip_string_and_ws('...'):
                if not self.skip_string_and_ws('('):
                    self.fail("Expecting '(' after 'sizeof...'.")
                if not self.match(identifier_re):
                    self.fail("Expecting identifier for 'sizeof...'.")
                ident = ASTIdentifier(self.matched_text)
                self.skip_ws()
                if not self.skip_string(")"):
                    self.fail("Expecting ')' to end 'sizeof...'.")
                return ASTSizeofParamPack(ident)
            if self.skip_string_and_ws('('):
                typ = self._parse_type(named=False)
                self.skip_ws()
                if not self.skip_string(')'):
                    self.fail("Expecting ')' to end 'sizeof'.")
                return ASTSizeofType(typ)
            expr = self._parse_unary_expression()
            return ASTSizeofExpr(expr)
        if self.skip_word_and_ws('alignof'):
            if not self.skip_string_and_ws('('):
                self.fail("Expecting '(' after 'alignof'.")
            typ = self._parse_type(named=False)
            self.skip_ws()
            if not self.skip_string(')'):
                self.fail("Expecting ')' to end 'alignof'.")
            return ASTAlignofExpr(typ)
        if self.skip_word_and_ws('noexcept'):
            if not self.skip_string_and_ws('('):
                self.fail("Expecting '(' after 'noexcept'.")
            expr = self._parse_expression()
            self.skip_ws()
            if not self.skip_string(')'):
                self.fail("Expecting ')' to end 'noexcept'.")
            return ASTNoexceptExpr(expr)
        # new-expression
        pos = self.pos
        rooted = self.skip_string('::')
        self.skip_ws()
        if not self.skip_word_and_ws('new'):
            self.pos = pos
        else:
            # new-placement[opt] new-type-id new-initializer[opt]
            # new-placement[opt] ( type-id ) new-initializer[opt]
            isNewTypeId = True
            if self.skip_string_and_ws('('):
                # either this is a new-placement or it's the second production
                # without placement, and it's actually the ( type-id ) part
                self.fail("Sorry, neither new-placement nor parenthesised type-id "
                          "in new-epression is supported yet.")
                # set isNewTypeId = False if it's (type-id)
            if isNewTypeId:
                declSpecs = self._parse_decl_specs(outer=None)
                decl = self._parse_declarator(named=False, paramMode="new")
            else:
                self.fail("Sorry, parenthesised type-id in new expression not yet supported.")
            lst = self._parse_expression_list_or_braced_init_list()
            return ASTNewExpr(rooted, isNewTypeId, ASTType(declSpecs, decl), lst)
        # delete-expression
        pos = self.pos
        rooted = self.skip_string('::')
        self.skip_ws()
        if not self.skip_word_and_ws('delete'):
            self.pos = pos
        else:
            array = self.skip_string_and_ws('[')
            if array and not self.skip_string_and_ws(']'):
                self.fail("Expected ']' in array delete-expression.")
            expr = self._parse_cast_expression()
            return ASTDeleteExpr(rooted, array, expr)
        return self._parse_postfix_expression()

    def _parse_cast_expression(self) -> ASTExpression:
        # -> unary  | "(" type-id ")" cast
        pos = self.pos
        self.skip_ws()
        if self.skip_string('('):
            try:
                typ = self._parse_type(False)
                if not self.skip_string(')'):
                    self.fail("Expected ')' in cast expression.")
                expr = self._parse_cast_expression()
                return ASTCastExpr(typ, expr)
            except DefinitionError as exCast:
                self.pos = pos
                try:
                    return self._parse_unary_expression()
                except DefinitionError as exUnary:
                    errs = []
                    errs.append((exCast, "If type cast expression"))
                    errs.append((exUnary, "If unary expression"))
                    raise self._make_multi_error(errs,
                                                 "Error in cast expression.") from exUnary
        else:
            return self._parse_unary_expression()

    def _parse_logical_or_expression(self, inTemplate: bool) -> ASTExpression:
        # logical-or     = logical-and      ||
        # logical-and    = inclusive-or     &&
        # inclusive-or   = exclusive-or     |
        # exclusive-or   = and              ^
        # and            = equality         &
        # equality       = relational       ==, !=
        # relational     = shift            <, >, <=, >=, <=>
        # shift          = additive         <<, >>
        # additive       = multiplicative   +, -
        # multiplicative = pm               *, /, %
        # pm             = cast             .*, ->*
        def _parse_bin_op_expr(self: DefinitionParser,
                               opId: int, inTemplate: bool) -> ASTExpression:
            if opId + 1 == len(_expression_bin_ops):
                def parser(inTemplate: bool) -> ASTExpression:
                    return self._parse_cast_expression()
            else:
                def parser(inTemplate: bool) -> ASTExpression:
                    return _parse_bin_op_expr(self, opId + 1, inTemplate=inTemplate)
            exprs = []
            ops = []
            exprs.append(parser(inTemplate=inTemplate))
            while True:
                self.skip_ws()
                if inTemplate and self.current_char == '>':
                    break
                pos = self.pos
                oneMore = False
                for op in _expression_bin_ops[opId]:
                    if op[0] in 'abcnox':
                        if not self.skip_word(op):
                            continue
                    else:
                        if not self.skip_string(op):
                            continue
                    if op == self.current_char == '&':
                        # don't split the && 'token'
                        self.pos -= 1
                        # and btw. && has lower precedence, so we are done
                        break
                    try:
                        expr = parser(inTemplate=inTemplate)
                        exprs.append(expr)
                        ops.append(op)
                        oneMore = True
                        break
                    except DefinitionError:
                        self.pos = pos
                if not oneMore:
                    break
            return ASTBinOpExpr(exprs, ops)
        return _parse_bin_op_expr(self, 0, inTemplate=inTemplate)

    def _parse_conditional_expression_tail(self, orExprHead: ASTExpression,
                                           inTemplate: bool) -> ASTConditionalExpr | None:
        # Consumes the orExprHead on success.

        # -> "?" expression ":" assignment-expression
        self.skip_ws()
        if not self.skip_string("?"):
            return None
        thenExpr = self._parse_expression()
        self.skip_ws()
        if not self.skip_string(":"):
            self.fail('Expected ":" after then-expression in conditional expression.')
        elseExpr = self._parse_assignment_expression(inTemplate)
        return ASTConditionalExpr(orExprHead, thenExpr, elseExpr)

    def _parse_assignment_expression(self, inTemplate: bool) -> ASTExpression:
        # -> conditional-expression
        #  | logical-or-expression assignment-operator initializer-clause
        #  | yield-expression -> "co_yield" assignment-expression
        #                      | "co_yield" braced-init-list
        #  | throw-expression -> "throw" assignment-expression[opt]
        # TODO: yield-expression
        # TODO: throw-expression

        # Now we have (after expanding conditional-expression:
        #     logical-or-expression
        #   | logical-or-expression "?" expression ":" assignment-expression
        #   | logical-or-expression assignment-operator initializer-clause
        leftExpr = self._parse_logical_or_expression(inTemplate=inTemplate)
        # the ternary operator
        condExpr = self._parse_conditional_expression_tail(leftExpr, inTemplate)
        if condExpr is not None:
            return condExpr
        # and actual assignment
        for op in _expression_assignment_ops:
            if op[0] in 'anox':
                if not self.skip_word(op):
                    continue
            else:
                if not self.skip_string(op):
                    continue
            rightExpr = self._parse_initializer_clause()
            return ASTAssignmentExpr(leftExpr, op, rightExpr)
        # just a logical-or-expression
        return leftExpr

    def _parse_constant_expression(self, inTemplate: bool) -> ASTExpression:
        # -> conditional-expression ->
        #    logical-or-expression
        #  | logical-or-expression "?" expression ":" assignment-expression
        orExpr = self._parse_logical_or_expression(inTemplate=inTemplate)
        condExpr = self._parse_conditional_expression_tail(orExpr, inTemplate)
        if condExpr is not None:
            return condExpr
        return orExpr

    def _parse_expression(self) -> ASTExpression:
        # -> assignment-expression
        #  | expression "," assignment-expression
        exprs = [self._parse_assignment_expression(inTemplate=False)]
        while True:
            self.skip_ws()
            if not self.skip_string(','):
                break
            exprs.append(self._parse_assignment_expression(inTemplate=False))
        if len(exprs) == 1:
            return exprs[0]
        else:
            return ASTCommaExpr(exprs)

    def _parse_expression_fallback(self, end: list[str],
                                   parser: Callable[[], ASTExpression],
                                   allow: bool = True) -> ASTExpression:
        # Stupidly "parse" an expression.
        # 'end' should be a list of characters which ends the expression.

        # first try to use the provided parser
        prevPos = self.pos
        try:
            return parser()
        except DefinitionError as e:
            # some places (e.g., template parameters) we really don't want to use fallback,
            # and for testing we may want to globally disable it
            if not allow or not self.allowFallbackExpressionParsing:
                raise
            self.warn("Parsing of expression failed. Using fallback parser."
                      " Error was:\n%s" % e)
            self.pos = prevPos
        # and then the fallback scanning
        assert end is not None
        self.skip_ws()
        startPos = self.pos
        if self.match(_string_re):
            value = self.matched_text
        else:
            # TODO: add handling of more bracket-like things, and quote handling
            brackets = {'(': ')', '{': '}', '[': ']', '<': '>'}
            symbols: list[str] = []
            while not self.eof:
                if (len(symbols) == 0 and self.current_char in end):
                    break
                if self.current_char in brackets:
                    symbols.append(brackets[self.current_char])
                elif len(symbols) > 0 and self.current_char == symbols[-1]:
                    symbols.pop()
                self.pos += 1
            if len(end) > 0 and self.eof:
                self.fail("Could not find end of expression starting at %d."
                          % startPos)
            value = self.definition[startPos:self.pos].strip()
        return ASTFallbackExpr(value.strip())

    # ==========================================================================

    def _parse_operator(self) -> ASTOperator:
        self.skip_ws()
        # adapted from the old code
        # yay, a regular operator definition
        if self.match(_operator_re):
            return ASTOperatorBuildIn(self.matched_text)

        # new/delete operator?
        for op in 'new', 'delete':
            if not self.skip_word(op):
                continue
            self.skip_ws()
            if self.skip_string('['):
                self.skip_ws()
                if not self.skip_string(']'):
                    self.fail('Expected "]" after  "operator ' + op + '["')
                op += '[]'
            return ASTOperatorBuildIn(op)

        # user-defined literal?
        if self.skip_string('""'):
            self.skip_ws()
            if not self.match(identifier_re):
                self.fail("Expected user-defined literal suffix.")
            identifier = ASTIdentifier(self.matched_text)
            return ASTOperatorLiteral(identifier)

        # oh well, looks like a cast operator definition.
        # In that case, eat another type.
        type = self._parse_type(named=False, outer="operatorCast")
        return ASTOperatorType(type)

    def _parse_template_argument_list(self) -> ASTTemplateArgs:
        # template-argument-list: (but we include the < and > here
        #    template-argument ...[opt]
        #    template-argument-list, template-argument ...[opt]
        # template-argument:
        #    constant-expression
        #    type-id
        #    id-expression
        self.skip_ws()
        if not self.skip_string_and_ws('<'):
            return None
        if self.skip_string('>'):
            return ASTTemplateArgs([], False)
        prevErrors = []
        templateArgs: list[ASTType | ASTTemplateArgConstant] = []
        packExpansion = False
        while 1:
            pos = self.pos
            parsedComma = False
            parsedEnd = False
            try:
                type = self._parse_type(named=False)
                self.skip_ws()
                if self.skip_string_and_ws('...'):
                    packExpansion = True
                    parsedEnd = True
                    if not self.skip_string('>'):
                        self.fail('Expected ">" after "..." in template argument list.')
                elif self.skip_string('>'):
                    parsedEnd = True
                elif self.skip_string(','):
                    parsedComma = True
                else:
                    self.fail('Expected "...>", ">" or "," in template argument list.')
                templateArgs.append(type)
            except DefinitionError as e:
                prevErrors.append((e, "If type argument"))
                self.pos = pos
                try:
                    value = self._parse_constant_expression(inTemplate=True)
                    self.skip_ws()
                    if self.skip_string_and_ws('...'):
                        packExpansion = True
                        parsedEnd = True
                        if not self.skip_string('>'):
                            self.fail('Expected ">" after "..." in template argument list.')
                    elif self.skip_string('>'):
                        parsedEnd = True
                    elif self.skip_string(','):
                        parsedComma = True
                    else:
                        self.fail('Expected "...>", ">" or "," in template argument list.')
                    templateArgs.append(ASTTemplateArgConstant(value))
                except DefinitionError as e:
                    self.pos = pos
                    prevErrors.append((e, "If non-type argument"))
                    header = "Error in parsing template argument list."
                    raise self._make_multi_error(prevErrors, header) from e
            if parsedEnd:
                assert not parsedComma
                break
            assert not packExpansion
        return ASTTemplateArgs(templateArgs, packExpansion)

    def _parse_nested_name(self, memberPointer: bool = False) -> ASTNestedName:
        names: list[ASTNestedNameElement] = []
        templates: list[bool] = []

        self.skip_ws()
        rooted = False
        if self.skip_string('::'):
            rooted = True
        while 1:
            self.skip_ws()
            if len(names) > 0:
                template = self.skip_word_and_ws('template')
            else:
                template = False
            templates.append(template)
            identOrOp: ASTIdentifier | ASTOperator | None = None
            if self.skip_word_and_ws('operator'):
                identOrOp = self._parse_operator()
            else:
                if not self.match(identifier_re):
                    if memberPointer and len(names) > 0:
                        templates.pop()
                        break
                    self.fail("Expected identifier in nested name.")
                identifier = self.matched_text
                # make sure there isn't a keyword
                if identifier in _keywords:
                    self.fail("Expected identifier in nested name, "
                              "got keyword: %s" % identifier)
                identOrOp = ASTIdentifier(identifier)
            # try greedily to get template arguments,
            # but otherwise a < might be because we are in an expression
            pos = self.pos
            try:
                templateArgs = self._parse_template_argument_list()
            except DefinitionError as ex:
                self.pos = pos
                templateArgs = None
                self.otherErrors.append(ex)
            names.append(ASTNestedNameElement(identOrOp, templateArgs))

            self.skip_ws()
            if not self.skip_string('::'):
                if memberPointer:
                    self.fail("Expected '::' in pointer to member (function).")
                break
        return ASTNestedName(names, templates, rooted)

    # ==========================================================================

    def _parse_simple_type_specifiers(self) -> ASTTrailingTypeSpecFundamental:
        modifier: str | None = None
        signedness: str | None = None
        width: list[str] = []
        typ: str | None = None
        names: list[str] = []  # the parsed sequence

        self.skip_ws()
        while self.match(_simple_type_specifiers_re):
            t = self.matched_text
            names.append(t)
            if t in ('auto', 'void', 'bool',
                     'char', 'wchar_t', 'char8_t', 'char16_t', 'char32_t',
                     'int', '__int64', '__int128',
                     'float', 'double',
                     '__float80', '_Float64x', '__float128', '_Float128'):
                if typ is not None:
                    self.fail(f"Can not have both {t} and {typ}.")
                typ = t
            elif t in ('signed', 'unsigned'):
                if signedness is not None:
                    self.fail(f"Can not have both {t} and {signedness}.")
                signedness = t
            elif t == 'short':
                if len(width) != 0:
                    self.fail(f"Can not have both {t} and {width[0]}.")
                width.append(t)
            elif t == 'long':
                if len(width) != 0 and width[0] != 'long':
                    self.fail(f"Can not have both {t} and {width[0]}.")
                width.append(t)
            elif t in ('_Imaginary', '_Complex'):
                if modifier is not None:
                    self.fail(f"Can not have both {t} and {modifier}.")
                modifier = t
            self.skip_ws()
        if len(names) == 0:
            return None

        if typ in ('auto', 'void', 'bool',
                   'wchar_t', 'char8_t', 'char16_t', 'char32_t',
                   '__float80', '_Float64x', '__float128', '_Float128'):
            if modifier is not None:
                self.fail(f"Can not have both {typ} and {modifier}.")
            if signedness is not None:
                self.fail(f"Can not have both {typ} and {signedness}.")
            if len(width) != 0:
                self.fail(f"Can not have both {typ} and {' '.join(width)}.")
        elif typ == 'char':
            if modifier is not None:
                self.fail(f"Can not have both {typ} and {modifier}.")
            if len(width) != 0:
                self.fail(f"Can not have both {typ} and {' '.join(width)}.")
        elif typ == 'int':
            if modifier is not None:
                self.fail(f"Can not have both {typ} and {modifier}.")
        elif typ in ('__int64', '__int128'):
            if modifier is not None:
                self.fail(f"Can not have both {typ} and {modifier}.")
            if len(width) != 0:
                self.fail(f"Can not have both {typ} and {' '.join(width)}.")
        elif typ == 'float':
            if signedness is not None:
                self.fail(f"Can not have both {typ} and {signedness}.")
            if len(width) != 0:
                self.fail(f"Can not have both {typ} and {' '.join(width)}.")
        elif typ == 'double':
            if signedness is not None:
                self.fail(f"Can not have both {typ} and {signedness}.")
            if len(width) > 1:
                self.fail(f"Can not have both {typ} and {' '.join(width)}.")
            if len(width) == 1 and width[0] != 'long':
                self.fail(f"Can not have both {typ} and {' '.join(width)}.")
        elif typ is None:
            if modifier is not None:
                self.fail(f"Can not have {modifier} without a floating point type.")
        else:
            msg = f'Unhandled type {typ}'
            raise AssertionError(msg)

        canonNames: list[str] = []
        if modifier is not None:
            canonNames.append(modifier)
        if signedness is not None:
            canonNames.append(signedness)
        canonNames.extend(width)
        if typ is not None:
            canonNames.append(typ)
        return ASTTrailingTypeSpecFundamental(names, canonNames)

    def _parse_trailing_type_spec(self) -> ASTTrailingTypeSpec:
        # fundamental types, https://en.cppreference.com/w/cpp/language/type
        # and extensions
        self.skip_ws()
        res = self._parse_simple_type_specifiers()
        if res is not None:
            return res

        # decltype
        self.skip_ws()
        if self.skip_word_and_ws('decltype'):
            if not self.skip_string_and_ws('('):
                self.fail("Expected '(' after 'decltype'.")
            if self.skip_word_and_ws('auto'):
                if not self.skip_string(')'):
                    self.fail("Expected ')' after 'decltype(auto'.")
                return ASTTrailingTypeSpecDecltypeAuto()
            expr = self._parse_expression()
            self.skip_ws()
            if not self.skip_string(')'):
                self.fail("Expected ')' after 'decltype(<expr>'.")
            return ASTTrailingTypeSpecDecltype(expr)

        # prefixed
        prefix = None
        self.skip_ws()
        for k in ('class', 'struct', 'enum', 'union', 'typename'):
            if self.skip_word_and_ws(k):
                prefix = k
                break
        nestedName = self._parse_nested_name()
        self.skip_ws()
        placeholderType = None
        if self.skip_word('auto'):
            placeholderType = 'auto'
        elif self.skip_word_and_ws('decltype'):
            if not self.skip_string_and_ws('('):
                self.fail("Expected '(' after 'decltype' in placeholder type specifier.")
            if not self.skip_word_and_ws('auto'):
                self.fail("Expected 'auto' after 'decltype(' in placeholder type specifier.")
            if not self.skip_string_and_ws(')'):
                self.fail("Expected ')' after 'decltype(auto' in placeholder type specifier.")
            placeholderType = 'decltype(auto)'
        return ASTTrailingTypeSpecName(prefix, nestedName, placeholderType)

    def _parse_parameters_and_qualifiers(
        self, paramMode: str,
    ) -> ASTParametersQualifiers | None:
        if paramMode == 'new':
            return None
        self.skip_ws()
        if not self.skip_string('('):
            if paramMode == 'function':
                self.fail('Expecting "(" in parameters-and-qualifiers.')
            else:
                return None
        args = []
        self.skip_ws()
        if not self.skip_string(')'):
            while 1:
                self.skip_ws()
                if self.skip_string('...'):
                    args.append(ASTFunctionParameter(None, True))
                    self.skip_ws()
                    if not self.skip_string(')'):
                        self.fail('Expected ")" after "..." in '
                                  'parameters-and-qualifiers.')
                    break
                # note: it seems that function arguments can always be named,
                # even in function pointers and similar.
                arg = self._parse_type_with_init(outer=None, named='single')
                # TODO: parse default parameters # TODO: didn't we just do that?
                args.append(ASTFunctionParameter(arg))

                self.skip_ws()
                if self.skip_string(','):
                    continue
                if self.skip_string(')'):
                    break
                self.fail('Expecting "," or ")" in parameters-and-qualifiers, '
                          f'got "{self.current_char}".')

        self.skip_ws()
        const = self.skip_word_and_ws('const')
        volatile = self.skip_word_and_ws('volatile')
        if not const:  # the can be permuted
            const = self.skip_word_and_ws('const')

        refQual = None
        if self.skip_string('&&'):
            refQual = '&&'
        if not refQual and self.skip_string('&'):
            refQual = '&'

        exceptionSpec = None
        self.skip_ws()
        if self.skip_string('noexcept'):
            if self.skip_string_and_ws('('):
                expr = self._parse_constant_expression(False)
                self.skip_ws()
                if not self.skip_string(')'):
                    self.fail("Expecting ')' to end 'noexcept'.")
                exceptionSpec = ASTNoexceptSpec(expr)
            else:
                exceptionSpec = ASTNoexceptSpec(None)

        self.skip_ws()
        if self.skip_string('->'):
            trailingReturn = self._parse_type(named=False)
        else:
            trailingReturn = None

        self.skip_ws()
        override = self.skip_word_and_ws('override')
        final = self.skip_word_and_ws('final')
        if not override:
            override = self.skip_word_and_ws(
                'override')  # they can be permuted

        attrs = self._parse_attribute_list()

        self.skip_ws()
        initializer = None
        # if this is a function pointer we should not swallow an initializer
        if paramMode == 'function' and self.skip_string('='):
            self.skip_ws()
            valid = ('0', 'delete', 'default')
            for w in valid:
                if self.skip_word_and_ws(w):
                    initializer = w
                    break
            if not initializer:
                self.fail(
                    'Expected "%s" in initializer-specifier.'
                    % '" or "'.join(valid))

        return ASTParametersQualifiers(
            args, volatile, const, refQual, exceptionSpec, trailingReturn,
            override, final, attrs, initializer)

    def _parse_decl_specs_simple(self, outer: str, typed: bool) -> ASTDeclSpecsSimple:
        """Just parse the simple ones."""
        storage = None
        threadLocal = None
        inline = None
        virtual = None
        explicitSpec = None
        consteval = None
        constexpr = None
        constinit = None
        volatile = None
        const = None
        friend = None
        attrs = []
        while 1:  # accept any permutation of a subset of some decl-specs
            self.skip_ws()
            if not const and typed:
                const = self.skip_word('const')
                if const:
                    continue
            if not volatile and typed:
                volatile = self.skip_word('volatile')
                if volatile:
                    continue
            if not storage:
                if outer in ('member', 'function'):
                    if self.skip_word('static'):
                        storage = 'static'
                        continue
                    if self.skip_word('extern'):
                        storage = 'extern'
                        continue
                if outer == 'member':
                    if self.skip_word('mutable'):
                        storage = 'mutable'
                        continue
                if self.skip_word('register'):
                    storage = 'register'
                    continue
            if not inline and outer in ('function', 'member'):
                inline = self.skip_word('inline')
                if inline:
                    continue
            if not constexpr and outer in ('member', 'function'):
                constexpr = self.skip_word("constexpr")
                if constexpr:
                    continue

            if outer == 'member':
                if not constinit:
                    constinit = self.skip_word('constinit')
                    if constinit:
                        continue
                if not threadLocal:
                    threadLocal = self.skip_word('thread_local')
                    if threadLocal:
                        continue
            if outer == 'function':
                if not consteval:
                    consteval = self.skip_word('consteval')
                    if consteval:
                        continue
                if not friend:
                    friend = self.skip_word('friend')
                    if friend:
                        continue
                if not virtual:
                    virtual = self.skip_word('virtual')
                    if virtual:
                        continue
                if not explicitSpec:
                    explicit = self.skip_word_and_ws('explicit')
                    if explicit:
                        expr: ASTExpression = None
                        if self.skip_string('('):
                            expr = self._parse_constant_expression(inTemplate=False)
                            if not expr:
                                self.fail("Expected constant expression after '('"
                                          " in explicit specifier.")
                            self.skip_ws()
                            if not self.skip_string(')'):
                                self.fail("Expected ')' to end explicit specifier.")
                        explicitSpec = ASTExplicitSpec(expr)
                        continue
            attr = self._parse_attribute()
            if attr:
                attrs.append(attr)
                continue
            break
        return ASTDeclSpecsSimple(storage, threadLocal, inline, virtual,
                                  explicitSpec, consteval, constexpr, constinit,
                                  volatile, const, friend, ASTAttributeList(attrs))

    def _parse_decl_specs(self, outer: str, typed: bool = True) -> ASTDeclSpecs:
        if outer:
            if outer not in ('type', 'member', 'function', 'templateParam'):
                raise Exception('Internal error, unknown outer "%s".' % outer)
        """
        storage-class-specifier function-specifier "constexpr"
        "volatile" "const" trailing-type-specifier

        storage-class-specifier ->
              "static" (only for member_object and function_object)
            | "register"

        function-specifier -> "inline" | "virtual" | "explicit" (only for
        function_object)

        "constexpr" (only for member_object and function_object)
        """
        leftSpecs = self._parse_decl_specs_simple(outer, typed)
        rightSpecs = None

        if typed:
            trailing = self._parse_trailing_type_spec()
            rightSpecs = self._parse_decl_specs_simple(outer, typed)
        else:
            trailing = None
        return ASTDeclSpecs(outer, leftSpecs, rightSpecs, trailing)

    def _parse_declarator_name_suffix(
        self, named: bool | str, paramMode: str, typed: bool,
    ) -> ASTDeclaratorNameParamQual | ASTDeclaratorNameBitField:
        # now we should parse the name, and then suffixes
        if named == 'maybe':
            pos = self.pos
            try:
                declId = self._parse_nested_name()
            except DefinitionError:
                self.pos = pos
                declId = None
        elif named == 'single':
            if self.match(identifier_re):
                identifier = ASTIdentifier(self.matched_text)
                nne = ASTNestedNameElement(identifier, None)
                declId = ASTNestedName([nne], [False], rooted=False)
                # if it's a member pointer, we may have '::', which should be an error
                self.skip_ws()
                if self.current_char == ':':
                    self.fail("Unexpected ':' after identifier.")
            else:
                declId = None
        elif named:
            declId = self._parse_nested_name()
        else:
            declId = None
        arrayOps = []
        while 1:
            self.skip_ws()
            if typed and self.skip_string('['):
                self.skip_ws()
                if self.skip_string(']'):
                    arrayOps.append(ASTArray(None))
                    continue

                def parser() -> ASTExpression:
                    return self._parse_expression()
                value = self._parse_expression_fallback([']'], parser)
                if not self.skip_string(']'):
                    self.fail("Expected ']' in end of array operator.")
                arrayOps.append(ASTArray(value))
                continue
            break
        paramQual = self._parse_parameters_and_qualifiers(paramMode)
        if paramQual is None and len(arrayOps) == 0:
            # perhaps a bit-field
            if named and paramMode == 'type' and typed:
                self.skip_ws()
                if self.skip_string(':'):
                    size = self._parse_constant_expression(inTemplate=False)
                    return ASTDeclaratorNameBitField(declId=declId, size=size)
        return ASTDeclaratorNameParamQual(declId=declId, arrayOps=arrayOps,
                                          paramQual=paramQual)

    def _parse_declarator(self, named: bool | str, paramMode: str,
                          typed: bool = True,
                          ) -> ASTDeclarator:
        # 'typed' here means 'parse return type stuff'
        if paramMode not in ('type', 'function', 'operatorCast', 'new'):
            raise Exception(
                "Internal error, unknown paramMode '%s'." % paramMode)
        prevErrors = []
        self.skip_ws()
        if typed and self.skip_string('*'):
            self.skip_ws()
            volatile = False
            const = False
            attrList = []
            while 1:
                if not volatile:
                    volatile = self.skip_word_and_ws('volatile')
                    if volatile:
                        continue
                if not const:
                    const = self.skip_word_and_ws('const')
                    if const:
                        continue
                attr = self._parse_attribute()
                if attr is not None:
                    attrList.append(attr)
                    continue
                break
            next = self._parse_declarator(named, paramMode, typed)
            return ASTDeclaratorPtr(next=next, volatile=volatile, const=const,
                                    attrs=ASTAttributeList(attrList))
        # TODO: shouldn't we parse an R-value ref here first?
        if typed and self.skip_string("&"):
            attrs = self._parse_attribute_list()
            next = self._parse_declarator(named, paramMode, typed)
            return ASTDeclaratorRef(next=next, attrs=attrs)
        if typed and self.skip_string("..."):
            next = self._parse_declarator(named, paramMode, False)
            return ASTDeclaratorParamPack(next=next)
        if typed and self.current_char == '(':  # note: peeking, not skipping
            if paramMode == "operatorCast":
                # TODO: we should be able to parse cast operators which return
                # function pointers. For now, just hax it and ignore.
                return ASTDeclaratorNameParamQual(declId=None, arrayOps=[],
                                                  paramQual=None)
            # maybe this is the beginning of params and quals,try that first,
            # otherwise assume it's noptr->declarator > ( ptr-declarator )
            pos = self.pos
            try:
                # assume this is params and quals
                res = self._parse_declarator_name_suffix(named, paramMode,
                                                         typed)
                return res
            except DefinitionError as exParamQual:
                prevErrors.append((exParamQual,
                                   "If declarator-id with parameters-and-qualifiers"))
                self.pos = pos
                try:
                    assert self.current_char == '('
                    self.skip_string('(')
                    # TODO: hmm, if there is a name, it must be in inner, right?
                    # TODO: hmm, if there must be parameters, they must be
                    #       inside, right?
                    inner = self._parse_declarator(named, paramMode, typed)
                    if not self.skip_string(')'):
                        self.fail("Expected ')' in \"( ptr-declarator )\"")
                    next = self._parse_declarator(named=False,
                                                  paramMode="type",
                                                  typed=typed)
                    return ASTDeclaratorParen(inner=inner, next=next)
                except DefinitionError as exNoPtrParen:
                    self.pos = pos
                    prevErrors.append((exNoPtrParen, "If parenthesis in noptr-declarator"))
                    header = "Error in declarator"
                    raise self._make_multi_error(prevErrors, header) from exNoPtrParen
        if typed:  # pointer to member
            pos = self.pos
            try:
                name = self._parse_nested_name(memberPointer=True)
                self.skip_ws()
                if not self.skip_string('*'):
                    self.fail("Expected '*' in pointer to member declarator.")
                self.skip_ws()
            except DefinitionError as e:
                self.pos = pos
                prevErrors.append((e, "If pointer to member declarator"))
            else:
                volatile = False
                const = False
                while 1:
                    if not volatile:
                        volatile = self.skip_word_and_ws('volatile')
                        if volatile:
                            continue
                    if not const:
                        const = self.skip_word_and_ws('const')
                        if const:
                            continue
                    break
                next = self._parse_declarator(named, paramMode, typed)
                return ASTDeclaratorMemPtr(name, const, volatile, next=next)
        pos = self.pos
        try:
            res = self._parse_declarator_name_suffix(named, paramMode, typed)
            # this is a heuristic for error messages, for when there is a < after a
            # nested name, but it was not a successful template argument list
            if self.current_char == '<':
                self.otherErrors.append(self._make_multi_error(prevErrors, ""))
            return res
        except DefinitionError as e:
            self.pos = pos
            prevErrors.append((e, "If declarator-id"))
            header = "Error in declarator or parameters-and-qualifiers"
            raise self._make_multi_error(prevErrors, header) from e

    def _parse_initializer(self, outer: str | None = None, allowFallback: bool = True,
                           ) -> ASTInitializer | None:
        # initializer                           # global vars
        # -> brace-or-equal-initializer
        #  | '(' expression-list ')'
        #
        # brace-or-equal-initializer            # member vars
        # -> '=' initializer-clause
        #  | braced-init-list
        #
        # initializer-clause  # function params, non-type template params (with '=' in front)
        # -> assignment-expression
        #  | braced-init-list
        #
        # we don't distinguish between global and member vars, so disallow paren:
        #
        # -> braced-init-list             # var only
        #  | '=' assignment-expression
        #  | '=' braced-init-list
        self.skip_ws()
        if outer == 'member':
            bracedInit = self._parse_braced_init_list()
            if bracedInit is not None:
                return ASTInitializer(bracedInit, hasAssign=False)

        if not self.skip_string('='):
            return None

        bracedInit = self._parse_braced_init_list()
        if bracedInit is not None:
            return ASTInitializer(bracedInit)

        if outer == 'member':
            fallbackEnd: list[str] = []
        elif outer == 'templateParam':
            fallbackEnd = [',', '>']
        elif outer is None:  # function parameter
            fallbackEnd = [',', ')']
        else:
            self.fail("Internal error, initializer for outer '%s' not "
                      "implemented." % outer)

        inTemplate = outer == 'templateParam'

        def parser() -> ASTExpression:
            return self._parse_assignment_expression(inTemplate=inTemplate)
        value = self._parse_expression_fallback(fallbackEnd, parser, allow=allowFallback)
        return ASTInitializer(value)

    def _parse_type(self, named: bool | str, outer: str | None = None) -> ASTType:
        """
        named=False|'maybe'|True: 'maybe' is e.g., for function objects which
        doesn't need to name the arguments

        outer == operatorCast: annoying case, we should not take the params
        """
        if outer:  # always named
            if outer not in ('type', 'member', 'function',
                             'operatorCast', 'templateParam'):
                raise Exception('Internal error, unknown outer "%s".' % outer)
            if outer != 'operatorCast':
                assert named
        if outer in ('type', 'function'):
            # We allow type objects to just be a name.
            # Some functions don't have normal return types: constructors,
            # destructors, cast operators
            prevErrors = []
            startPos = self.pos
            # first try without the type
            try:
                declSpecs = self._parse_decl_specs(outer=outer, typed=False)
                decl = self._parse_declarator(named=True, paramMode=outer,
                                              typed=False)
                mustEnd = True
                if outer == 'function':
                    # Allow trailing requires on functions.
                    self.skip_ws()
                    if re.compile(r'requires\b').match(self.definition, self.pos):
                        mustEnd = False
                if mustEnd:
                    self.assert_end(allowSemicolon=True)
            except DefinitionError as exUntyped:
                if outer == 'type':
                    desc = "If just a name"
                elif outer == 'function':
                    desc = "If the function has no return type"
                else:
                    raise AssertionError from exUntyped
                prevErrors.append((exUntyped, desc))
                self.pos = startPos
                try:
                    declSpecs = self._parse_decl_specs(outer=outer)
                    decl = self._parse_declarator(named=True, paramMode=outer)
                except DefinitionError as exTyped:
                    self.pos = startPos
                    if outer == 'type':
                        desc = "If typedef-like declaration"
                    elif outer == 'function':
                        desc = "If the function has a return type"
                    else:
                        raise AssertionError from exUntyped
                    prevErrors.append((exTyped, desc))
                    # Retain the else branch for easier debugging.
                    # TODO: it would be nice to save the previous stacktrace
                    #       and output it here.
                    if True:
                        if outer == 'type':
                            header = "Type must be either just a name or a "
                            header += "typedef-like declaration."
                        elif outer == 'function':
                            header = "Error when parsing function declaration."
                        else:
                            raise AssertionError from exUntyped
                        raise self._make_multi_error(prevErrors, header) from exTyped
                    else:  # NoQA: RET506
                        # For testing purposes.
                        # do it again to get the proper traceback (how do you
                        # reliably save a traceback when an exception is
                        # constructed?)
                        self.pos = startPos
                        typed = True
                        declSpecs = self._parse_decl_specs(outer=outer, typed=typed)
                        decl = self._parse_declarator(named=True, paramMode=outer,
                                                      typed=typed)
        else:
            paramMode = 'type'
            if outer == 'member':
                named = True
            elif outer == 'operatorCast':
                paramMode = 'operatorCast'
                outer = None
            elif outer == 'templateParam':
                named = 'single'
            declSpecs = self._parse_decl_specs(outer=outer)
            decl = self._parse_declarator(named=named, paramMode=paramMode)
        return ASTType(declSpecs, decl)

    def _parse_type_with_init(
            self, named: bool | str,
            outer: str) -> ASTTypeWithInit | ASTTemplateParamConstrainedTypeWithInit:
        if outer:
            assert outer in ('type', 'member', 'function', 'templateParam')
        type = self._parse_type(outer=outer, named=named)
        if outer != 'templateParam':
            init = self._parse_initializer(outer=outer)
            return ASTTypeWithInit(type, init)
        # it could also be a constrained type parameter, e.g., C T = int&
        pos = self.pos
        eExpr = None
        try:
            init = self._parse_initializer(outer=outer, allowFallback=False)
            # note: init may be None if there is no =
            if init is None:
                return ASTTypeWithInit(type, None)
            # we parsed an expression, so we must have a , or a >,
            # otherwise the expression didn't get everything
            self.skip_ws()
            if self.current_char != ',' and self.current_char != '>':
                # pretend it didn't happen
                self.pos = pos
                init = None
            else:
                # we assume that it was indeed an expression
                return ASTTypeWithInit(type, init)
        except DefinitionError as e:
            self.pos = pos
            eExpr = e
        if not self.skip_string("="):
            return ASTTypeWithInit(type, None)
        try:
            typeInit = self._parse_type(named=False, outer=None)
            return ASTTemplateParamConstrainedTypeWithInit(type, typeInit)
        except DefinitionError as eType:
            if eExpr is None:
                raise
            errs = []
            errs.append((eExpr, "If default template argument is an expression"))
            errs.append((eType, "If default template argument is a type"))
            msg = "Error in non-type template parameter"
            msg += " or constrained template parameter."
            raise self._make_multi_error(errs, msg) from eType

    def _parse_type_using(self) -> ASTTypeUsing:
        name = self._parse_nested_name()
        self.skip_ws()
        if not self.skip_string('='):
            return ASTTypeUsing(name, None)
        type = self._parse_type(False, None)
        return ASTTypeUsing(name, type)

    def _parse_concept(self) -> ASTConcept:
        nestedName = self._parse_nested_name()
        self.skip_ws()
        initializer = self._parse_initializer('member')
        return ASTConcept(nestedName, initializer)

    def _parse_class(self) -> ASTClass:
        attrs = self._parse_attribute_list()
        name = self._parse_nested_name()
        self.skip_ws()
        final = self.skip_word_and_ws('final')
        bases = []
        self.skip_ws()
        if self.skip_string(':'):
            while 1:
                self.skip_ws()
                visibility = None
                virtual = False
                pack = False
                if self.skip_word_and_ws('virtual'):
                    virtual = True
                if self.match(_visibility_re):
                    visibility = self.matched_text
                    self.skip_ws()
                if not virtual and self.skip_word_and_ws('virtual'):
                    virtual = True
                baseName = self._parse_nested_name()
                self.skip_ws()
                pack = self.skip_string('...')
                bases.append(ASTBaseClass(baseName, visibility, virtual, pack))
                self.skip_ws()
                if self.skip_string(','):
                    continue
                break
        return ASTClass(name, final, bases, attrs)

    def _parse_union(self) -> ASTUnion:
        attrs = self._parse_attribute_list()
        name = self._parse_nested_name()
        return ASTUnion(name, attrs)

    def _parse_enum(self) -> ASTEnum:
        scoped = None  # is set by CPPEnumObject
        attrs = self._parse_attribute_list()
        name = self._parse_nested_name()
        self.skip_ws()
        underlyingType = None
        if self.skip_string(':'):
            underlyingType = self._parse_type(named=False)
        return ASTEnum(name, scoped, underlyingType, attrs)

    def _parse_enumerator(self) -> ASTEnumerator:
        name = self._parse_nested_name()
        attrs = self._parse_attribute_list()
        self.skip_ws()
        init = None
        if self.skip_string('='):
            self.skip_ws()

            def parser() -> ASTExpression:
                return self._parse_constant_expression(inTemplate=False)
            initVal = self._parse_expression_fallback([], parser)
            init = ASTInitializer(initVal)
        return ASTEnumerator(name, init, attrs)

    # ==========================================================================

    def _parse_template_parameter(self) -> ASTTemplateParam:
        self.skip_ws()
        if self.skip_word('template'):
            # declare a template template parameter
            nestedParams = self._parse_template_parameter_list()
        else:
            nestedParams = None

        pos = self.pos
        try:
            # Unconstrained type parameter or template type parameter
            key = None
            self.skip_ws()
            if self.skip_word_and_ws('typename'):
                key = 'typename'
            elif self.skip_word_and_ws('class'):
                key = 'class'
            elif nestedParams:
                self.fail("Expected 'typename' or 'class' after "
                          "template template parameter list.")
            else:
                self.fail("Expected 'typename' or 'class' in the "
                          "beginning of template type parameter.")
            self.skip_ws()
            parameterPack = self.skip_string('...')
            self.skip_ws()
            if self.match(identifier_re):
                identifier = ASTIdentifier(self.matched_text)
            else:
                identifier = None
            self.skip_ws()
            if not parameterPack and self.skip_string('='):
                default = self._parse_type(named=False, outer=None)
            else:
                default = None
                if self.current_char not in ',>':
                    self.fail('Expected "," or ">" after (template) type parameter.')
            data = ASTTemplateKeyParamPackIdDefault(key, identifier,
                                                    parameterPack, default)
            if nestedParams:
                return ASTTemplateParamTemplateType(nestedParams, data)
            else:
                return ASTTemplateParamType(data)
        except DefinitionError as eType:
            if nestedParams:
                raise
            try:
                # non-type parameter or constrained type parameter
                self.pos = pos
                param = self._parse_type_with_init('maybe', 'templateParam')
                self.skip_ws()
                parameterPack = self.skip_string('...')
                return ASTTemplateParamNonType(param, parameterPack)
            except DefinitionError as eNonType:
                self.pos = pos
                header = "Error when parsing template parameter."
                errs = []
                errs.append(
                    (eType, "If unconstrained type parameter or template type parameter"))
                errs.append(
                    (eNonType, "If constrained type parameter or non-type parameter"))
                raise self._make_multi_error(errs, header) from None

    def _parse_template_parameter_list(self) -> ASTTemplateParams:
        # only: '<' parameter-list '>'
        # we assume that 'template' has just been parsed
        templateParams: list[ASTTemplateParam] = []
        self.skip_ws()
        if not self.skip_string("<"):
            self.fail("Expected '<' after 'template'")
        while 1:
            pos = self.pos
            err = None
            try:
                param = self._parse_template_parameter()
                templateParams.append(param)
            except DefinitionError as eParam:
                self.pos = pos
                err = eParam
            self.skip_ws()
            if self.skip_string('>'):
                requiresClause = self._parse_requires_clause()
                return ASTTemplateParams(templateParams, requiresClause)
            elif self.skip_string(','):
                continue
            else:
                header = "Error in template parameter list."
                errs = []
                if err:
                    errs.append((err, "If parameter"))
                try:
                    self.fail('Expected "," or ">".')
                except DefinitionError as e:
                    errs.append((e, "If no parameter"))
                logger.debug(errs)
                raise self._make_multi_error(errs, header)

    def _parse_template_introduction(self) -> ASTTemplateIntroduction | None:
        pos = self.pos
        try:
            concept = self._parse_nested_name()
        except Exception:
            self.pos = pos
            return None
        self.skip_ws()
        if not self.skip_string('{'):
            self.pos = pos
            return None

        # for sure it must be a template introduction now
        params = []
        while 1:
            self.skip_ws()
            parameterPack = self.skip_string('...')
            self.skip_ws()
            if not self.match(identifier_re):
                self.fail("Expected identifier in template introduction list.")
            txt_identifier = self.matched_text
            # make sure there isn't a keyword
            if txt_identifier in _keywords:
                self.fail("Expected identifier in template introduction list, "
                          "got keyword: %s" % txt_identifier)
            identifier = ASTIdentifier(txt_identifier)
            params.append(ASTTemplateIntroductionParameter(identifier, parameterPack))

            self.skip_ws()
            if self.skip_string('}'):
                break
            if self.skip_string(','):
                continue
            self.fail('Error in template introduction list. Expected ",", or "}".')
        return ASTTemplateIntroduction(concept, params)

    def _parse_requires_clause(self) -> ASTRequiresClause | None:
        # requires-clause -> 'requires' constraint-logical-or-expression
        # constraint-logical-or-expression
        #   -> constraint-logical-and-expression
        #    | constraint-logical-or-expression '||' constraint-logical-and-expression
        # constraint-logical-and-expression
        #   -> primary-expression
        #    | constraint-logical-and-expression '&&' primary-expression
        self.skip_ws()
        if not self.skip_word('requires'):
            return None

        def parse_and_expr(self: DefinitionParser) -> ASTExpression:
            andExprs = []
            ops = []
            andExprs.append(self._parse_primary_expression())
            while True:
                self.skip_ws()
                oneMore = False
                if self.skip_string('&&'):
                    oneMore = True
                    ops.append('&&')
                elif self.skip_word('and'):
                    oneMore = True
                    ops.append('and')
                if not oneMore:
                    break
                andExprs.append(self._parse_primary_expression())
            if len(andExprs) == 1:
                return andExprs[0]
            else:
                return ASTBinOpExpr(andExprs, ops)

        orExprs = []
        ops = []
        orExprs.append(parse_and_expr(self))
        while True:
            self.skip_ws()
            oneMore = False
            if self.skip_string('||'):
                oneMore = True
                ops.append('||')
            elif self.skip_word('or'):
                oneMore = True
                ops.append('or')
            if not oneMore:
                break
            orExprs.append(parse_and_expr(self))
        if len(orExprs) == 1:
            return ASTRequiresClause(orExprs[0])
        else:
            return ASTRequiresClause(ASTBinOpExpr(orExprs, ops))

    def _parse_template_declaration_prefix(self, objectType: str,
                                           ) -> ASTTemplateDeclarationPrefix | None:
        templates: list[ASTTemplateParams | ASTTemplateIntroduction] = []
        while 1:
            self.skip_ws()
            # the saved position is only used to provide a better error message
            params: ASTTemplateParams | ASTTemplateIntroduction | None = None
            pos = self.pos
            if self.skip_word("template"):
                try:
                    params = self._parse_template_parameter_list()
                except DefinitionError as e:
                    if objectType == 'member' and len(templates) == 0:
                        return ASTTemplateDeclarationPrefix(None)
                    else:
                        raise e
                if objectType == 'concept' and params.requiresClause is not None:
                    self.fail('requires-clause not allowed for concept')
            else:
                params = self._parse_template_introduction()
                if not params:
                    break
            if objectType == 'concept' and len(templates) > 0:
                self.pos = pos
                self.fail("More than 1 template parameter list for concept.")
            templates.append(params)
        if len(templates) == 0 and objectType == 'concept':
            self.fail('Missing template parameter list for concept.')
        if len(templates) == 0:
            return None
        else:
            return ASTTemplateDeclarationPrefix(templates)

    def _check_template_consistency(self, nestedName: ASTNestedName,
                                    templatePrefix: ASTTemplateDeclarationPrefix,
                                    fullSpecShorthand: bool, isMember: bool = False,
                                    ) -> ASTTemplateDeclarationPrefix:
        numArgs = nestedName.num_templates()
        isMemberInstantiation = False
        if not templatePrefix:
            numParams = 0
        else:
            if isMember and templatePrefix.templates is None:
                numParams = 0
                isMemberInstantiation = True
            else:
                numParams = len(templatePrefix.templates)
        if numArgs + 1 < numParams:
            self.fail("Too few template argument lists compared to parameter"
                      " lists. Argument lists: %d, Parameter lists: %d."
                      % (numArgs, numParams))
        if numArgs > numParams:
            numExtra = numArgs - numParams
            if not fullSpecShorthand and not isMemberInstantiation:
                msg = (
                    f'Too many template argument lists compared to parameter lists. '
                    f'Argument lists: {numArgs:d}, Parameter lists: {numParams:d}, '
                    f'Extra empty parameters lists prepended: {numExtra:d}. '
                    'Declaration:\n\t'
                )
                if templatePrefix:
                    msg += f"{templatePrefix}\n\t"
                msg += str(nestedName)
                self.warn(msg)

            newTemplates: list[ASTTemplateParams | ASTTemplateIntroduction] = [
                ASTTemplateParams([], requiresClause=None)
                for _i in range(numExtra)
            ]
            if templatePrefix and not isMemberInstantiation:
                newTemplates.extend(templatePrefix.templates)
            templatePrefix = ASTTemplateDeclarationPrefix(newTemplates)
        return templatePrefix

    def parse_declaration(self, objectType: str, directiveType: str) -> ASTDeclaration:
        if objectType not in ('class', 'union', 'function', 'member', 'type',
                              'concept', 'enum', 'enumerator'):
            raise Exception('Internal error, unknown objectType "%s".' % objectType)
        if directiveType not in ('class', 'struct', 'union', 'function', 'member', 'var',
                                 'type', 'concept',
                                 'enum', 'enum-struct', 'enum-class', 'enumerator'):
            raise Exception('Internal error, unknown directiveType "%s".' % directiveType)
        visibility = None
        templatePrefix = None
        trailingRequiresClause = None
        declaration: Any = None

        self.skip_ws()
        if self.match(_visibility_re):
            visibility = self.matched_text

        if objectType in ('type', 'concept', 'member', 'function', 'class', 'union'):
            templatePrefix = self._parse_template_declaration_prefix(objectType)

        if objectType == 'type':
            prevErrors = []
            pos = self.pos
            try:
                if not templatePrefix:
                    declaration = self._parse_type(named=True, outer='type')
            except DefinitionError as e:
                prevErrors.append((e, "If typedef-like declaration"))
                self.pos = pos
            pos = self.pos
            try:
                if not declaration:
                    declaration = self._parse_type_using()
            except DefinitionError as e:
                self.pos = pos
                prevErrors.append((e, "If type alias or template alias"))
                header = "Error in type declaration."
                raise self._make_multi_error(prevErrors, header) from e
        elif objectType == 'concept':
            declaration = self._parse_concept()
        elif objectType == 'member':
            declaration = self._parse_type_with_init(named=True, outer='member')
        elif objectType == 'function':
            declaration = self._parse_type(named=True, outer='function')
            trailingRequiresClause = self._parse_requires_clause()
        elif objectType == 'class':
            declaration = self._parse_class()
        elif objectType == 'union':
            declaration = self._parse_union()
        elif objectType == 'enum':
            declaration = self._parse_enum()
        elif objectType == 'enumerator':
            declaration = self._parse_enumerator()
        else:
            raise AssertionError
        templatePrefix = self._check_template_consistency(declaration.name,
                                                          templatePrefix,
                                                          fullSpecShorthand=False,
                                                          isMember=objectType == 'member')
        self.skip_ws()
        semicolon = self.skip_string(';')
        return ASTDeclaration(objectType, directiveType, visibility,
                              templatePrefix, declaration,
                              trailingRequiresClause, semicolon)

    def parse_namespace_object(self) -> ASTNamespace:
        templatePrefix = self._parse_template_declaration_prefix(objectType="namespace")
        name = self._parse_nested_name()
        templatePrefix = self._check_template_consistency(name, templatePrefix,
                                                          fullSpecShorthand=False)
        res = ASTNamespace(name, templatePrefix)
        res.objectType = 'namespace'  # type: ignore[attr-defined]
        return res

    def parse_xref_object(self) -> tuple[ASTNamespace | ASTDeclaration, bool]:
        pos = self.pos
        try:
            templatePrefix = self._parse_template_declaration_prefix(objectType="xref")
            name = self._parse_nested_name()
            # if there are '()' left, just skip them
            self.skip_ws()
            self.skip_string('()')
            self.assert_end()
            templatePrefix = self._check_template_consistency(name, templatePrefix,
                                                              fullSpecShorthand=True)
            res1 = ASTNamespace(name, templatePrefix)
            res1.objectType = 'xref'  # type: ignore[attr-defined]
            return res1, True
        except DefinitionError as e1:
            try:
                self.pos = pos
                res2 = self.parse_declaration('function', 'function')
                # if there are '()' left, just skip them
                self.skip_ws()
                self.skip_string('()')
                self.assert_end()
                return res2, False
            except DefinitionError as e2:
                errs = []
                errs.append((e1, "If shorthand ref"))
                errs.append((e2, "If full function ref"))
                msg = "Error in cross-reference."
                raise self._make_multi_error(errs, msg) from e2

    def parse_expression(self) -> ASTExpression | ASTType:
        pos = self.pos
        try:
            expr = self._parse_expression()
            self.skip_ws()
            self.assert_end()
            return expr
        except DefinitionError as exExpr:
            self.pos = pos
            try:
                typ = self._parse_type(False)
                self.skip_ws()
                self.assert_end()
                return typ
            except DefinitionError as exType:
                header = "Error when parsing (type) expression."
                errs = []
                errs.append((exExpr, "If expression"))
                errs.append((exType, "If type"))
                raise self._make_multi_error(errs, header) from exType
