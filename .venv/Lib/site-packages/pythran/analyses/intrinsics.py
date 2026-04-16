""" Intrinsics gathers all intrinsics referenced in a module. """

from pythran.passmanager import ModuleAnalysis
import pythran.intrinsic as intrinsic
from pythran.utils import attr_to_path


class Intrinsics(ModuleAnalysis):

    """ Gather all intrinsics used in the module
    """

    """ Result is a set of intrinsic values. """
    ResultType = set

    def visit_Attribute(self, node):
        obj, _ = attr_to_path(node)
        if isinstance(obj, intrinsic.Intrinsic):
            self.result.add(obj)
        self.generic_visit(node)
