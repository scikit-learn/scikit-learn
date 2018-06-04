from rope.base.oi.type_hinting.providers import interfaces


class ParamProvider(interfaces.IParamProvider):

    def __init__(self, *delegates):
        """
        :type delegates: list[rope.base.oi.type_hinting.providers.interfaces.IParamProvider]
        """
        self._delegates = delegates

    def __call__(self, pyfunc, param_name):
        """
        :type pyfunc: rope.base.pyobjectsdef.PyFunction
        :type param_name: str
        :rtype: rope.base.pyobjects.PyDefinedObject | rope.base.pyobjects.PyObject or None
        """
        for delegate in self._delegates:
            result = delegate(pyfunc, param_name)
            if result:
                return result


class ReturnProvider(interfaces.IReturnProvider):

    def __init__(self, *delegates):
        """
        :type delegates: list[rope.base.oi.type_hinting.providers.interfaces.IReturnProvider]
        """
        self._delegates = delegates

    def __call__(self, pyfunc):
        """
        :type pyfunc: rope.base.pyobjectsdef.PyFunction
        :rtype: rope.base.pyobjects.PyDefinedObject | rope.base.pyobjects.PyObject or None
        """
        for delegate in self._delegates:
            result = delegate(pyfunc)
            if result:
                return result


class AssignmentProvider(interfaces.IAssignmentProvider):

    def __init__(self, *delegates):
        """
        :type delegates: list[rope.base.oi.type_hinting.providers.interfaces.IAssignmentProvider]
        """
        self._delegates = delegates

    def __call__(self, pyname):
        """
        :type pyname: rope.base.pynamesdef.AssignedName
        :rtype: rope.base.pyobjects.PyDefinedObject | rope.base.pyobjects.PyObject or None
        """
        for delegate in self._delegates:
            result = delegate(pyname)
            if result:
                return result
