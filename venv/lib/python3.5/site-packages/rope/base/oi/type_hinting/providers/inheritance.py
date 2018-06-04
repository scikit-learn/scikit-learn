from rope.base.oi.type_hinting import utils
from rope.base.oi.type_hinting.providers import interfaces


class ParamProvider(interfaces.IParamProvider):

    def __init__(self, delegate):
        """
        :type delegate: rope.base.oi.type_hinting.providers.interfaces.IParamProvider
        """
        self._delegate = delegate

    def __call__(self, pyfunc, param_name):
        """
        :type pyfunc: rope.base.pyobjectsdef.PyFunction
        :type param_name: str
        :rtype: rope.base.pyobjects.PyDefinedObject | rope.base.pyobjects.PyObject or None
        """
        superfunc = pyfunc
        while superfunc:
            result = self._delegate(superfunc, param_name)
            if result:
                return result
            superfunc = utils.get_super_func(superfunc)


class ReturnProvider(interfaces.IReturnProvider):

    def __init__(self, delegate):
        """
        :type delegate: rope.base.oi.type_hinting.providers.interfaces.IReturnProvider
        """
        self._delegate = delegate

    def __call__(self, pyfunc):
        """
        :type pyfunc: rope.base.pyobjectsdef.PyFunction
        :rtype: rope.base.pyobjects.PyDefinedObject | rope.base.pyobjects.PyObject or None
        """
        superfunc = pyfunc
        while superfunc:
            result = self._delegate(superfunc)
            if result:
                return result
            superfunc = utils.get_super_func(superfunc)


class AssignmentProvider(interfaces.IAssignmentProvider):

    def __init__(self, delegate):
        """
        :type delegate: rope.base.oi.type_hinting.providers.interfaces.IAssignmentProvider
        """
        self._delegate = delegate

    def __call__(self, pyname):
        """
        :type pyname: rope.base.pynamesdef.AssignedName
        :rtype: rope.base.pyobjects.PyDefinedObject | rope.base.pyobjects.PyObject or None
        """
        super_pyname = pyname
        while super_pyname:
            result = self._delegate(super_pyname)
            if result:
                return result
            super_pyname = utils.get_super_assignment(super_pyname)
