from inspect import signature, Parameter
from functools import wraps
import warnings


# This function was adapted from scikit-learn
# github.com/scikit-learn/scikit-learn/blob/master/sklearn/utils/validation.py
def _deprecate_positional_args(f):
    """Decorator for methods that issues warnings for positional arguments.

    Using the keyword-only argument syntax in pep 3102, arguments after the
    * will issue a warning when passed as a positional argument.

    Parameters
    ----------
    f : function
        function to check arguments on

    """
    sig = signature(f)
    kwonly_args = []
    all_args = []

    for name, param in sig.parameters.items():
        if param.kind == Parameter.POSITIONAL_OR_KEYWORD:
            all_args.append(name)
        elif param.kind == Parameter.KEYWORD_ONLY:
            kwonly_args.append(name)

    @wraps(f)
    def inner_f(*args, **kwargs):
        extra_args = len(args) - len(all_args)
        if extra_args > 0:
            plural = "s" if extra_args > 1 else ""
            article = "" if plural else "a "
            warnings.warn(
                "Pass the following variable{} as {}keyword arg{}: {}. "
                "From version 0.12, the only valid positional argument "
                "will be `data`, and passing other arguments without an "
                "explicit keyword will result in an error or misinterpretation."
                .format(plural, article, plural,
                        ", ".join(kwonly_args[:extra_args])),
                FutureWarning
            )
        kwargs.update({k: arg for k, arg in zip(sig.parameters, args)})
        return f(**kwargs)
    return inner_f


def share_init_params_with_map(cls):
    """Make cls.map a classmethod with same signature as cls.__init__."""
    map_sig = signature(cls.map)
    init_sig = signature(cls.__init__)

    new = [v for k, v in init_sig.parameters.items() if k != "self"]
    new.insert(0, map_sig.parameters["cls"])
    cls.map.__signature__ = map_sig.replace(parameters=new)
    cls.map.__doc__ = cls.__init__.__doc__

    cls.map = classmethod(cls.map)

    return cls
