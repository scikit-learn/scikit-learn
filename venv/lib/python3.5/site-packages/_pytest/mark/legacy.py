"""
this is a place where we put datastructures used by legacy apis
we hope ot remove
"""
import attr
import keyword

from . import MarkInfo, MarkDecorator

from _pytest.config import UsageError


@attr.s
class MarkMapping(object):
    """Provides a local mapping for markers where item access
    resolves to True if the marker is present. """

    own_mark_names = attr.ib()

    @classmethod
    def from_keywords(cls, keywords):
        mark_names = set()
        for key, value in keywords.items():
            if isinstance(value, MarkInfo) or isinstance(value, MarkDecorator):
                mark_names.add(key)
        return cls(mark_names)

    def __getitem__(self, name):
        return name in self.own_mark_names


class KeywordMapping(object):
    """Provides a local mapping for keywords.
    Given a list of names, map any substring of one of these names to True.
    """

    def __init__(self, names):
        self._names = names

    @classmethod
    def from_item(cls, item):
        mapped_names = set()

        # Add the names of the current item and any parent items
        import pytest
        for item in item.listchain():
            if not isinstance(item, pytest.Instance):
                mapped_names.add(item.name)

        # Add the names added as extra keywords to current or parent items
        for name in item.listextrakeywords():
            mapped_names.add(name)

        # Add the names attached to the current function through direct assignment
        if hasattr(item, 'function'):
            for name in item.function.__dict__:
                mapped_names.add(name)

        return cls(mapped_names)

    def __getitem__(self, subname):
        for name in self._names:
            if subname in name:
                return True
        return False


python_keywords_allowed_list = ["or", "and", "not"]


def matchmark(colitem, markexpr):
    """Tries to match on any marker names, attached to the given colitem."""
    return eval(markexpr, {}, MarkMapping.from_keywords(colitem.keywords))


def matchkeyword(colitem, keywordexpr):
    """Tries to match given keyword expression to given collector item.

    Will match on the name of colitem, including the names of its parents.
    Only matches names of items which are either a :class:`Class` or a
    :class:`Function`.
    Additionally, matches on names in the 'extra_keyword_matches' set of
    any item, as well as names directly assigned to test functions.
    """
    mapping = KeywordMapping.from_item(colitem)
    if " " not in keywordexpr:
        # special case to allow for simple "-k pass" and "-k 1.3"
        return mapping[keywordexpr]
    elif keywordexpr.startswith("not ") and " " not in keywordexpr[4:]:
        return not mapping[keywordexpr[4:]]
    for kwd in keywordexpr.split():
        if keyword.iskeyword(kwd) and kwd not in python_keywords_allowed_list:
            raise UsageError("Python keyword '{}' not accepted in expressions passed to '-k'".format(kwd))
    try:
        return eval(keywordexpr, {}, mapping)
    except SyntaxError:
        raise UsageError("Wrong expression passed to '-k': {}".format(keywordexpr))
