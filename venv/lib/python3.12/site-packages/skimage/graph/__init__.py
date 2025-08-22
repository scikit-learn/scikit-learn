"""
Graph-based operations, e.g., shortest paths.

This includes creating adjacency graphs of pixels in an image, finding the
central pixel in an image, finding (minimum-cost) paths across pixels, merging
and cutting of graphs, etc.

"""

import lazy_loader as _lazy

__getattr__, __dir__, __all__ = _lazy.attach_stub(__name__, __file__)
