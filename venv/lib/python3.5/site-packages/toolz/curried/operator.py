from __future__ import absolute_import

import operator

from toolz.functoolz import curry, num_required_args, has_keywords


def should_curry(f):
    num = num_required_args(f)
    return num is None or num > 1 or num == 1 and has_keywords(f) is not False


locals().update(
    dict((name, curry(f) if should_curry(f) else f)
         for name, f in vars(operator).items() if callable(f)),
)

# Clean up the namespace.
del curry
del num_required_args
del has_keywords
del operator
del should_curry
