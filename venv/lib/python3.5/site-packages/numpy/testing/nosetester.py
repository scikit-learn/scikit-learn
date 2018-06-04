"""
Back compatibility nosetester module. It will import the appropriate
set of tools

"""
from .nose_tools.nosetester import *

__all__ = ['get_package_name', 'run_module_suite', 'NoseTester',
           '_numpy_tester', 'get_package_name', 'import_nose',
           'suppress_warnings']
