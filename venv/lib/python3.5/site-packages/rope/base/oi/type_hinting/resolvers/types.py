from rope.base.oi.type_hinting import evaluate
from rope.base.oi.type_hinting.resolvers import interfaces


class Resolver(interfaces.IResolver):
    def __call__(self, hint, pyobject):
        """
        :param hint: For example "List[int]" or "(Foo, Bar) -> Baz" or simple "Foo"
        :type hint: str
        :type pyobject: rope.base.pyobjects.PyDefinedObject | rope.base.pyobjects.PyObject
        :rtype: rope.base.pyobjects.PyDefinedObject | rope.base.pyobjects.PyObject or None
        """
        try:
            return evaluate.evaluate(hint, pyobject)
        except (Exception):
            pass
