"""
****
YAML
****
Read and write NetworkX graphs in YAML format.

"YAML is a data serialization format designed for human readability 
and interaction with scripting languages."
See http://www.yaml.org for documentation.

Format
------
http://pyyaml.org/wiki/PyYAML

"""
__author__ = """Aric Hagberg (hagberg@lanl.gov)"""
#    Copyright (C) 2004-2018 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.

__all__ = ['read_yaml', 'write_yaml']

import networkx as nx
from networkx.utils import open_file


@open_file(1, mode='w')
def write_yaml(G_to_be_yaml, path_for_yaml_output, **kwds):
    """Write graph G in YAML format to path. 

    YAML is a data serialization format designed for human readability 
    and interaction with scripting languages [1]_.

    Parameters
    ----------
    G : graph
       A NetworkX graph
    path : file or string
       File or filename to write. 
       Filenames ending in .gz or .bz2 will be compressed.

    Notes
    -----
    To use encoding on the output file include e.g. `encoding='utf-8'`
    in the keyword arguments.

    Examples
    --------
    >>> G=nx.path_graph(4)
    >>> nx.write_yaml(G,'test.yaml')

    References
    ----------
    .. [1] http://www.yaml.org
    """
    try:
        import yaml
    except ImportError:
        raise ImportError("write_yaml() requires PyYAML: http://pyyaml.org/")
    yaml.dump(G_to_be_yaml, path_for_yaml_output, **kwds)


@open_file(0, mode='r')
def read_yaml(path):
    """Read graph in YAML format from path.

    YAML is a data serialization format designed for human readability 
    and interaction with scripting languages [1]_.

    Parameters
    ----------
    path : file or string
       File or filename to read.  Filenames ending in .gz or .bz2 
       will be uncompressed.

    Returns
    -------
    G : NetworkX graph

    Examples
    --------
    >>> G=nx.path_graph(4)
    >>> nx.write_yaml(G,'test.yaml')
    >>> G=nx.read_yaml('test.yaml')

    References
    ----------
    .. [1] http://www.yaml.org

    """
    try:
        import yaml
    except ImportError:
        raise ImportError("read_yaml() requires PyYAML: http://pyyaml.org/")

    G = yaml.load(path)
    return G


# fixture for nose tests
def setup_module(module):
    from nose import SkipTest
    try:
        import yaml
    except:
        raise SkipTest("PyYAML not available")

# fixture for nose tests


def teardown_module(module):
    import os
    os.unlink('test.yaml')
