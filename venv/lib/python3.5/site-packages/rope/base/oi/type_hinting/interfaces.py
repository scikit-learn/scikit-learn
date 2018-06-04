class ITypeHintingFactory(object):

    def make_param_provider(self):
        """
        :rtype: rope.base.oi.type_hinting.providers.interfaces.IParamProvider
        """
        raise NotImplementedError

    def make_return_provider(self):
        """
        :rtype: rope.base.oi.type_hinting.providers.interfaces.IReturnProvider
        """
        raise NotImplementedError

    def make_assignment_provider(self):
        """
        :rtype: rope.base.oi.type_hinting.providers.interfaces.IAssignmentProvider
        """
        raise NotImplementedError

    def make_resolver(self):
        """
        :rtype: rope.base.oi.type_hinting.resolvers.interfaces.IResolver
        """
        raise NotImplementedError
