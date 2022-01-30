import warnings
from seaborn import *  # noqa
reset_orig()  # noqa

msg = (
    "As seaborn no longer sets a default style on import, the seaborn.apionly "
    "module is deprecated. It will be removed in a future version."
)
warnings.warn(msg, UserWarning)
