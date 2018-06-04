import pydoc

from jedi.evaluate.utils import ignored
from jedi.evaluate.filters import AbstractNameDefinition

try:
    from pydoc_data import topics as pydoc_topics
except ImportError:
    # Python 2
    try:
        import pydoc_topics
    except ImportError:
        # This is for Python 3 embeddable version, which dont have
        # pydoc_data module in its file python3x.zip.
        pydoc_topics = None


def get_operator(evaluator, string, pos):
    return Keyword(evaluator, string, pos)


class KeywordName(AbstractNameDefinition):
    api_type = u'keyword'

    def __init__(self, evaluator, name):
        self.evaluator = evaluator
        self.string_name = name
        self.parent_context = evaluator.builtins_module

    def infer(self):
        return [Keyword(self.evaluator, self.string_name, (0, 0))]


class Keyword(object):
    api_type = u'keyword'

    def __init__(self, evaluator, name, pos):
        self.name = KeywordName(evaluator, name)
        self.start_pos = pos
        self.parent = evaluator.builtins_module

    @property
    def names(self):
        """ For a `parsing.Name` like comparision """
        return [self.name]

    def py__doc__(self, include_call_signature=False):
        return imitate_pydoc(self.name.string_name)

    def __repr__(self):
        return '<%s: %s>' % (type(self).__name__, self.name)


def imitate_pydoc(string):
    """
    It's not possible to get the pydoc's without starting the annoying pager
    stuff.
    """
    if pydoc_topics is None:
        return ''

    # str needed because of possible unicode stuff in py2k (pydoc doesn't work
    # with unicode strings)
    string = str(string)
    h = pydoc.help
    with ignored(KeyError):
        # try to access symbols
        string = h.symbols[string]
        string, _, related = string.partition(' ')

    get_target = lambda s: h.topics.get(s, h.keywords.get(s))
    while isinstance(string, str):
        string = get_target(string)

    try:
        # is a tuple now
        label, related = string
    except TypeError:
        return ''

    try:
        return pydoc_topics.topics[label].strip() if pydoc_topics else ''
    except KeyError:
        return ''
