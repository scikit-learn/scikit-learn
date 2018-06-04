from warnings import warn

warn("IPython.utils.pickleutil has moved to ipykernel.pickleutil", stacklevel=2)

from ipykernel.pickleutil import *
