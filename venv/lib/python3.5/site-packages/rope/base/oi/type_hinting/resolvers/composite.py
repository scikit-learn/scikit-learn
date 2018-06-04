from rope.base.oi.type_hinting.resolvers import interfaces


class Resolver(interfaces.IResolver):

    def __init__(self, *delegates):
        """
        :type delegates: list[rope.base.oi.type_hinting.resolvers.interfaces.IResolver]
        """
        self._delegates = delegates

    def __call__(self, hint, pyobject):
        """
        :param hint: For example "List[int]" or "(Foo, Bar) -> Baz" or simple "Foo"
        :type hint: str
        :type pyobject: rope.base.pyobjects.PyDefinedObject | rope.base.pyobjects.PyObject
        :rtype: rope.base.pyobjects.PyDefinedObject | rope.base.pyobjects.PyObject or None
        """
        for delegate in self._delegates:
            result = delegate(hint, pyobject)
            if result:
                return result
