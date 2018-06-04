from __future__ import absolute_import, division, print_function

from collections import deque

from dask.core import istask, subs


def head(task):
    """Return the top level node of a task"""

    if istask(task):
        return task[0]
    elif isinstance(task, list):
        return list
    else:
        return task


def args(task):
    """Get the arguments for the current task"""

    if istask(task):
        return task[1:]
    elif isinstance(task, list):
        return task
    else:
        return ()


class Traverser(object):
    """Traverser interface for tasks.

    Class for storing the state while performing a preorder-traversal of a
    task.

    Parameters
    ----------
    term : task
        The task to be traversed

    Attributes
    ----------
    term
        The current element in the traversal
    current
        The head of the current element in the traversal. This is simply `head`
        applied to the attribute `term`.
    """

    def __init__(self, term, stack=None):
        self.term = term
        if not stack:
            self._stack = deque([END])
        else:
            self._stack = stack

    def __iter__(self):
        while self.current is not END:
            yield self.current
            self.next()

    def copy(self):
        """Copy the traverser in its current state.

        This allows the traversal to be pushed onto a stack, for easy
        backtracking."""

        return Traverser(self.term, deque(self._stack))

    def next(self):
        """Proceed to the next term in the preorder traversal."""

        subterms = args(self.term)
        if not subterms:
            # No subterms, pop off stack
            self.term = self._stack.pop()
        else:
            self.term = subterms[0]
            self._stack.extend(reversed(subterms[1:]))

    @property
    def current(self):
        return head(self.term)

    def skip(self):
        """Skip over all subterms of the current level in the traversal"""
        self.term = self._stack.pop()


class Token(object):
    """A token object.

    Used to express certain objects in the traversal of a task or pattern."""

    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return self.name


# A variable to represent *all* variables in a discrimination net
VAR = Token('?')
# Represents the end of the traversal of an expression. We can't use `None`,
# 'False', etc... here, as anything may be an argument to a function.
END = Token('end')


class Node(tuple):
    """A Discrimination Net node."""

    __slots__ = ()

    def __new__(cls, edges=None, patterns=None):
        edges = edges if edges else {}
        patterns = patterns if patterns else []
        return tuple.__new__(cls, (edges, patterns))

    @property
    def edges(self):
        """A dictionary, where the keys are edges, and the values are nodes"""
        return self[0]

    @property
    def patterns(self):
        """A list of all patterns that currently match at this node"""
        return self[1]


class RewriteRule(object):
    """A rewrite rule.

    Expresses `lhs` -> `rhs`, for variables `vars`.

    Parameters
    ----------
    lhs : task
        The left-hand-side of the rewrite rule.
    rhs : task or function
        The right-hand-side of the rewrite rule. If it's a task, variables in
        `rhs` will be replaced by terms in the subject that match the variables
        in `lhs`. If it's a function, the function will be called with a dict
        of such matches.
    vars: tuple, optional
        Tuple of variables found in the lhs. Variables can be represented as
        any hashable object; a good convention is to use strings. If there are
        no variables, this can be omitted.

    Examples
    --------
    Here's a `RewriteRule` to replace all nested calls to `list`, so that
    `(list, (list, 'x'))` is replaced with `(list, 'x')`, where `'x'` is a
    variable.

    >>> lhs = (list, (list, 'x'))
    >>> rhs = (list, 'x')
    >>> variables = ('x',)
    >>> rule = RewriteRule(lhs, rhs, variables)

    Here's a more complicated rule that uses a callable right-hand-side. A
    callable `rhs` takes in a dictionary mapping variables to their matching
    values. This rule replaces all occurrences of `(list, 'x')` with `'x'` if
    `'x'` is a list itself.

    >>> lhs = (list, 'x')
    >>> def repl_list(sd):
    ...     x = sd['x']
    ...     if isinstance(x, list):
    ...         return x
    ...     else:
    ...         return (list, x)
    >>> rule = RewriteRule(lhs, repl_list, variables)
    """

    def __init__(self, lhs, rhs, vars=()):
        if not isinstance(vars, tuple):
            raise TypeError("vars must be a tuple of variables")
        self.lhs = lhs
        if callable(rhs):
            self.subs = rhs
        else:
            self.subs = self._apply
        self.rhs = rhs
        self._varlist = [t for t in Traverser(lhs) if t in vars]
        # Reduce vars down to just variables found in lhs
        self.vars = tuple(sorted(set(self._varlist)))

    def _apply(self, sub_dict):
        term = self.rhs
        for key, val in sub_dict.items():
            term = subs(term, key, val)
        return term

    def __str__(self):
        return "RewriteRule({0}, {1}, {2})".format(self.lhs, self.rhs,
                                                   self.vars)

    def __repr__(self):
        return str(self)


class RuleSet(object):
    """A set of rewrite rules.

    Forms a structure for fast rewriting over a set of rewrite rules. This
    allows for syntactic matching of terms to patterns for many patterns at
    the same time.

    Examples
    --------

    >>> def f(*args): pass
    >>> def g(*args): pass
    >>> def h(*args): pass
    >>> from operator import add

    >>> rs = RuleSet(                 # Make RuleSet with two Rules
    ...         RewriteRule((add, 'x', 0), 'x', ('x',)),
    ...         RewriteRule((f, (g, 'x'), 'y'),
    ...                     (h, 'x', 'y'),
    ...                     ('x', 'y')))

    >>> rs.rewrite((add, 2, 0))       # Apply ruleset to single task
    2

    >>> rs.rewrite((f, (g, 'a', 3)))  # doctest: +SKIP
    (h, 'a', 3)

    >>> dsk = {'a': (add, 2, 0),      # Apply ruleset to full dask graph
    ...        'b': (f, (g, 'a', 3))}

    >>> from toolz import valmap
    >>> valmap(rs.rewrite, dsk)  # doctest: +SKIP
    {'a': 2,
     'b': (h, 'a', 3)}

    Attributes
    ----------
    rules : list
        A list of `RewriteRule`s included in the `RuleSet`.
    """

    def __init__(self, *rules):
        """Create a `RuleSet` for a number of rules

        Parameters
        ----------
        rules
            One or more instances of RewriteRule
        """
        self._net = Node()
        self.rules = []
        for p in rules:
            self.add(p)

    def add(self, rule):
        """Add a rule to the RuleSet.

        Parameters
        ----------
        rule : RewriteRule
        """

        if not isinstance(rule, RewriteRule):
            raise TypeError("rule must be instance of RewriteRule")
        vars = rule.vars
        curr_node = self._net
        ind = len(self.rules)
        # List of variables, in order they appear in the POT of the term
        for t in Traverser(rule.lhs):
            prev_node = curr_node
            if t in vars:
                t = VAR
            if t in curr_node.edges:
                curr_node = curr_node.edges[t]
            else:
                curr_node.edges[t] = Node()
                curr_node = curr_node.edges[t]
        # We've reached a leaf node. Add the term index to this leaf.
        prev_node.edges[t].patterns.append(ind)
        self.rules.append(rule)

    def iter_matches(self, term):
        """A generator that lazily finds matchings for term from the RuleSet.

        Parameters
        ----------
        term : task

        Yields
        ------
        Tuples of `(rule, subs)`, where `rule` is the rewrite rule being
        matched, and `subs` is a dictionary mapping the variables in the lhs
        of the rule to their matching values in the term."""

        S = Traverser(term)
        for m, syms in _match(S, self._net):
            for i in m:
                rule = self.rules[i]
                subs = _process_match(rule, syms)
                if subs is not None:
                    yield rule, subs

    def _rewrite(self, term):
        """Apply the rewrite rules in RuleSet to top level of term"""

        for rule, sd in self.iter_matches(term):
            # We use for (...) because it's fast in all cases for getting the
            # first element from the match iterator. As we only want that
            # element, we break here
            term = rule.subs(sd)
            break
        return term

    def rewrite(self, task, strategy="bottom_up"):
        """Apply the `RuleSet` to `task`.

        This applies the most specific matching rule in the RuleSet to the
        task, using the provided strategy.

        Parameters
        ----------
        term: a task
            The task to be rewritten
        strategy: str, optional
            The rewriting strategy to use. Options are "bottom_up" (default),
            or "top_level".

        Examples
        --------
        Suppose there was a function `add` that returned the sum of 2 numbers,
        and another function `double` that returned twice its input:

        >>> add = lambda x, y: x + y
        >>> double = lambda x: 2*x

        Now suppose `double` was *significantly* faster than `add`, so
        you'd like to replace all expressions `(add, x, x)` with `(double,
        x)`, where `x` is a variable. This can be expressed as a rewrite rule:

        >>> rule = RewriteRule((add, 'x', 'x'), (double, 'x'), ('x',))
        >>> rs = RuleSet(rule)

        This can then be applied to terms to perform the rewriting:

        >>> term = (add, (add, 2, 2), (add, 2, 2))
        >>> rs.rewrite(term)  # doctest: +SKIP
        (double, (double, 2))

        If we only wanted to apply this to the top level of the term, the
        `strategy` kwarg can be set to "top_level".

        >>> rs.rewrite(term)  # doctest: +SKIP
        (double, (add, 2, 2))
        """
        return strategies[strategy](self, task)


def _top_level(net, term):
    return net._rewrite(term)


def _bottom_up(net, term):
    if istask(term):
        term = (head(term),) + tuple(_bottom_up(net, t) for t in args(term))
    elif isinstance(term, list):
        term = [_bottom_up(net, t) for t in args(term)]
    return net._rewrite(term)


strategies = {'top_level': _top_level,
              'bottom_up': _bottom_up}


def _match(S, N):
    """Structural matching of term S to discrimination net node N."""

    stack = deque()
    restore_state_flag = False
    # matches are stored in a tuple, because all mutations result in a copy,
    # preventing operations from changing matches stored on the stack.
    matches = ()
    while True:
        if S.current is END:
            yield N.patterns, matches
        try:
            # This try-except block is to catch hashing errors from un-hashable
            # types. This allows for variables to be matched with un-hashable
            # objects.
            n = N.edges.get(S.current, None)
            if n and not restore_state_flag:
                stack.append((S.copy(), N, matches))
                N = n
                S.next()
                continue
        except TypeError:
            pass
        n = N.edges.get(VAR, None)
        if n:
            restore_state_flag = False
            matches = matches + (S.term,)
            S.skip()
            N = n
            continue
        try:
            # Backtrack here
            (S, N, matches) = stack.pop()
            restore_state_flag = True
        except Exception:
            return


def _process_match(rule, syms):
    """Process a match to determine if it is correct, and to find the correct
    substitution that will convert the term into the pattern.

    Parameters
    ----------
    rule : RewriteRule
    syms : iterable
        Iterable of subterms that match a corresponding variable.

    Returns
    -------
    A dictionary of {vars : subterms} describing the substitution to make the
    pattern equivalent with the term. Returns `None` if the match is
    invalid."""

    subs = {}
    varlist = rule._varlist
    if not len(varlist) == len(syms):
        raise RuntimeError("length of varlist doesn't match length of syms.")
    for v, s in zip(varlist, syms):
        if v in subs and subs[v] != s:
            return None
        else:
            subs[v] = s
    return subs
