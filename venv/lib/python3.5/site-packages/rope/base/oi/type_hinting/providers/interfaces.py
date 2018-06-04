class IParamProvider(object):

    def __call__(self, pyfunc, param_name):
        """
        :type pyfunc: rope.base.pyobjectsdef.PyFunction
        :type param_name: str
        :rtype: rope.base.pyobjects.PyDefinedObject | rope.base.pyobjects.PyObject or None
        """
        raise NotImplementedError


class IReturnProvider(object):
    """
    :type resolve: rope.base.oi.type_hinting.resolvers.interfaces.IResolver
    """
    resolve = None

    def __call__(self, pyfunc):
        """
        :type pyfunc: rope.base.pyobjectsdef.PyFunction
        :rtype: rope.base.pyobjects.PyDefinedObject | rope.base.pyobjects.PyObject or None
        """
        raise NotImplementedError


class IAssignmentProvider(object):
    """
    :type resolve: rope.base.oi.type_hinting.resolvers.interfaces.IResolver
    """
    resolve = None

    def __call__(self, pyname):
        """
        :type pyname: rope.base.pynamesdef.AssignedName
        :rtype: rope.base.pyobjects.PyDefinedObject | rope.base.pyobjects.PyObject or None
        """
        raise NotImplementedError
