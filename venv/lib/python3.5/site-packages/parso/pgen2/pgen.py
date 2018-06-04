# Copyright 2004-2005 Elemental Security, Inc. All Rights Reserved.
# Licensed to PSF under a Contributor Agreement.

# Modifications:
# Copyright 2014 David Halter. Integration into Jedi.
# Modifications are dual-licensed: MIT and PSF.

from parso.pgen2 import grammar
from parso.python import token
from parso.python import tokenize
from parso.utils import parse_version_string


class ParserGenerator(object):
    def __init__(self, bnf_text, token_namespace):
        self._bnf_text = bnf_text
        self.generator = tokenize.tokenize(
            bnf_text,
            version_info=parse_version_string('3.6')
        )
        self._gettoken()  # Initialize lookahead
        self.dfas, self.startsymbol = self._parse()
        self.first = {}  # map from symbol name to set of tokens
        self._addfirstsets()
        self._token_namespace = token_namespace

    def make_grammar(self):
        c = grammar.Grammar(self._bnf_text)
        names = list(self.dfas.keys())
        names.sort()
        # TODO do we still need this?
        names.remove(self.startsymbol)
        names.insert(0, self.startsymbol)
        for name in names:
            i = 256 + len(c.symbol2number)
            c.symbol2number[name] = i
            c.number2symbol[i] = name
        for name in names:
            dfa = self.dfas[name]
            states = []
            for state in dfa:
                arcs = []
                for label, next in state.arcs.items():
                    arcs.append((self._make_label(c, label), dfa.index(next)))
                if state.isfinal:
                    arcs.append((0, dfa.index(state)))
                states.append(arcs)
            c.states.append(states)
            c.dfas[c.symbol2number[name]] = (states, self._make_first(c, name))
        c.start = c.symbol2number[self.startsymbol]
        return c

    def _make_first(self, c, name):
        rawfirst = self.first[name]
        first = {}
        for label in rawfirst:
            ilabel = self._make_label(c, label)
            ##assert ilabel not in first # XXX failed on <> ... !=
            first[ilabel] = 1
        return first

    def _make_label(self, c, label):
        # XXX Maybe this should be a method on a subclass of converter?
        ilabel = len(c.labels)
        if label[0].isalpha():
            # Either a symbol name or a named token
            if label in c.symbol2number:
                # A symbol name (a non-terminal)
                if label in c.symbol2label:
                    return c.symbol2label[label]
                else:
                    c.labels.append((c.symbol2number[label], None))
                    c.symbol2label[label] = ilabel
                    c.label2symbol[ilabel] = label
                    return ilabel
            else:
                # A named token (NAME, NUMBER, STRING)
                itoken = getattr(self._token_namespace, label, None)
                assert isinstance(itoken, int), label
                if itoken in c.tokens:
                    return c.tokens[itoken]
                else:
                    c.labels.append((itoken, None))
                    c.tokens[itoken] = ilabel
                    return ilabel
        else:
            # Either a keyword or an operator
            assert label[0] in ('"', "'"), label
            value = eval(label)
            if value[0].isalpha():
                # A keyword
                if value in c.keywords:
                    return c.keywords[value]
                else:
                    # TODO this might be an issue?! Using token.NAME here?
                    c.labels.append((token.NAME, value))
                    c.keywords[value] = ilabel
                    return ilabel
            else:
                # An operator (any non-numeric token)
                itoken = self._token_namespace.generate_token_id(value)
                if itoken in c.tokens:
                    return c.tokens[itoken]
                else:
                    c.labels.append((itoken, None))
                    c.tokens[itoken] = ilabel
                    return ilabel

    def _addfirstsets(self):
        names = list(self.dfas.keys())
        names.sort()
        for name in names:
            if name not in self.first:
                self._calcfirst(name)
            #print name, self.first[name].keys()

    def _calcfirst(self, name):
        dfa = self.dfas[name]
        self.first[name] = None  # dummy to detect left recursion
        state = dfa[0]
        totalset = {}
        overlapcheck = {}
        for label, next in state.arcs.items():
            if label in self.dfas:
                if label in self.first:
                    fset = self.first[label]
                    if fset is None:
                        raise ValueError("recursion for rule %r" % name)
                else:
                    self._calcfirst(label)
                    fset = self.first[label]
                totalset.update(fset)
                overlapcheck[label] = fset
            else:
                totalset[label] = 1
                overlapcheck[label] = {label: 1}
        inverse = {}
        for label, itsfirst in overlapcheck.items():
            for symbol in itsfirst:
                if symbol in inverse:
                    raise ValueError("rule %s is ambiguous; %s is in the"
                                     " first sets of %s as well as %s" %
                                     (name, symbol, label, inverse[symbol]))
                inverse[symbol] = label
        self.first[name] = totalset

    def _parse(self):
        dfas = {}
        startsymbol = None
        # MSTART: (NEWLINE | RULE)* ENDMARKER
        while self.type != token.ENDMARKER:
            while self.type == token.NEWLINE:
                self._gettoken()
            # RULE: NAME ':' RHS NEWLINE
            name = self._expect(token.NAME)
            self._expect(token.COLON)
            a, z = self._parse_rhs()
            self._expect(token.NEWLINE)
            #self._dump_nfa(name, a, z)
            dfa = self._make_dfa(a, z)
            #self._dump_dfa(name, dfa)
            # oldlen = len(dfa)
            self._simplify_dfa(dfa)
            # newlen = len(dfa)
            dfas[name] = dfa
            #print name, oldlen, newlen
            if startsymbol is None:
                startsymbol = name
        return dfas, startsymbol

    def _make_dfa(self, start, finish):
        # To turn an NFA into a DFA, we define the states of the DFA
        # to correspond to *sets* of states of the NFA.  Then do some
        # state reduction.  Let's represent sets as dicts with 1 for
        # values.
        assert isinstance(start, NFAState)
        assert isinstance(finish, NFAState)

        def closure(state):
            base = {}
            addclosure(state, base)
            return base

        def addclosure(state, base):
            assert isinstance(state, NFAState)
            if state in base:
                return
            base[state] = 1
            for label, next in state.arcs:
                if label is None:
                    addclosure(next, base)

        states = [DFAState(closure(start), finish)]
        for state in states:  # NB states grows while we're iterating
            arcs = {}
            for nfastate in state.nfaset:
                for label, next in nfastate.arcs:
                    if label is not None:
                        addclosure(next, arcs.setdefault(label, {}))
            for label, nfaset in arcs.items():
                for st in states:
                    if st.nfaset == nfaset:
                        break
                else:
                    st = DFAState(nfaset, finish)
                    states.append(st)
                state.addarc(st, label)
        return states  # List of DFAState instances; first one is start

    def _dump_nfa(self, name, start, finish):
        print("Dump of NFA for", name)
        todo = [start]
        for i, state in enumerate(todo):
            print("  State", i, state is finish and "(final)" or "")
            for label, next in state.arcs:
                if next in todo:
                    j = todo.index(next)
                else:
                    j = len(todo)
                    todo.append(next)
                if label is None:
                    print("    -> %d" % j)
                else:
                    print("    %s -> %d" % (label, j))

    def _dump_dfa(self, name, dfa):
        print("Dump of DFA for", name)
        for i, state in enumerate(dfa):
            print("  State", i, state.isfinal and "(final)" or "")
            for label, next in state.arcs.items():
                print("    %s -> %d" % (label, dfa.index(next)))

    def _simplify_dfa(self, dfa):
        # This is not theoretically optimal, but works well enough.
        # Algorithm: repeatedly look for two states that have the same
        # set of arcs (same labels pointing to the same nodes) and
        # unify them, until things stop changing.

        # dfa is a list of DFAState instances
        changes = True
        while changes:
            changes = False
            for i, state_i in enumerate(dfa):
                for j in range(i + 1, len(dfa)):
                    state_j = dfa[j]
                    if state_i == state_j:
                        #print "  unify", i, j
                        del dfa[j]
                        for state in dfa:
                            state.unifystate(state_j, state_i)
                        changes = True
                        break

    def _parse_rhs(self):
        # RHS: ALT ('|' ALT)*
        a, z = self._parse_alt()
        if self.value != "|":
            return a, z
        else:
            aa = NFAState()
            zz = NFAState()
            aa.addarc(a)
            z.addarc(zz)
            while self.value == "|":
                self._gettoken()
                a, z = self._parse_alt()
                aa.addarc(a)
                z.addarc(zz)
            return aa, zz

    def _parse_alt(self):
        # ALT: ITEM+
        a, b = self._parse_item()
        while (self.value in ("(", "[") or
               self.type in (token.NAME, token.STRING)):
            c, d = self._parse_item()
            b.addarc(c)
            b = d
        return a, b

    def _parse_item(self):
        # ITEM: '[' RHS ']' | ATOM ['+' | '*']
        if self.value == "[":
            self._gettoken()
            a, z = self._parse_rhs()
            self._expect(token.RSQB)
            a.addarc(z)
            return a, z
        else:
            a, z = self._parse_atom()
            value = self.value
            if value not in ("+", "*"):
                return a, z
            self._gettoken()
            z.addarc(a)
            if value == "+":
                return a, z
            else:
                return a, a

    def _parse_atom(self):
        # ATOM: '(' RHS ')' | NAME | STRING
        if self.value == "(":
            self._gettoken()
            a, z = self._parse_rhs()
            self._expect(token.RPAR)
            return a, z
        elif self.type in (token.NAME, token.STRING):
            a = NFAState()
            z = NFAState()
            a.addarc(z, self.value)
            self._gettoken()
            return a, z
        else:
            self._raise_error("expected (...) or NAME or STRING, got %s/%s",
                              self.type, self.value)

    def _expect(self, type):
        if self.type != type:
            self._raise_error("expected %s(%s), got %s(%s)",
                              type, token.tok_name[type], self.type, self.value)
        value = self.value
        self._gettoken()
        return value

    def _gettoken(self):
        tup = next(self.generator)
        while tup[0] in (token.COMMENT, token.NL):
            tup = next(self.generator)
        self.type, self.value, self.begin, prefix = tup

    def _raise_error(self, msg, *args):
        if args:
            try:
                msg = msg % args
            except:
                msg = " ".join([msg] + list(map(str, args)))
        line = self._bnf_text.splitlines()[self.begin[0] - 1]
        raise SyntaxError(msg, ('<grammar>', self.begin[0],
                                self.begin[1], line))


class NFAState(object):
    def __init__(self):
        self.arcs = []  # list of (label, NFAState) pairs

    def addarc(self, next, label=None):
        assert label is None or isinstance(label, str)
        assert isinstance(next, NFAState)
        self.arcs.append((label, next))


class DFAState(object):
    def __init__(self, nfaset, final):
        assert isinstance(nfaset, dict)
        assert isinstance(next(iter(nfaset)), NFAState)
        assert isinstance(final, NFAState)
        self.nfaset = nfaset
        self.isfinal = final in nfaset
        self.arcs = {}  # map from label to DFAState

    def addarc(self, next, label):
        assert isinstance(label, str)
        assert label not in self.arcs
        assert isinstance(next, DFAState)
        self.arcs[label] = next

    def unifystate(self, old, new):
        for label, next in self.arcs.items():
            if next is old:
                self.arcs[label] = new

    def __eq__(self, other):
        # Equality test -- ignore the nfaset instance variable
        assert isinstance(other, DFAState)
        if self.isfinal != other.isfinal:
            return False
        # Can't just return self.arcs == other.arcs, because that
        # would invoke this method recursively, with cycles...
        if len(self.arcs) != len(other.arcs):
            return False
        for label, next in self.arcs.items():
            if next is not other.arcs.get(label):
                return False
        return True

    __hash__ = None  # For Py3 compatibility.


def generate_grammar(bnf_text, token_namespace):
    """
    ``bnf_text`` is a grammar in extended BNF (using * for repetition, + for
    at-least-once repetition, [] for optional parts, | for alternatives and ()
    for grouping).

    It's not EBNF according to ISO/IEC 14977. It's a dialect Python uses in its
    own parser.
    """
    p = ParserGenerator(bnf_text, token_namespace)
    return p.make_grammar()
