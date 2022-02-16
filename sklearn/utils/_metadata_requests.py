"""
Metadata Routing Utility
"""

# Author: Adrin Jalali <adrin.jalali@gmail.com>
# License: BSD 3 clause

import inspect
from copy import deepcopy
from enum import Enum, EnumMeta
from warnings import warn
from collections import namedtuple
from typing import Union, Optional
from ._bunch import Bunch

# This namedtuple is used to store a (mapping, routing) pair. Mapping is a
# MethodMapping object, and routing is the output of `get_metadata_routing`.
# MetadataRouter stores a collection of these namedtuples.
RouterMappingPair = namedtuple("RouterMappingPair", ["mapping", "router"])

# A namedtuple storing a single method route. A collection of these namedtuples
# is stored in a MetadataRouter.
MethodPair = namedtuple("MethodPair", ["callee", "caller"])


class MemberCheckEnumMeta(EnumMeta):
    """Enum metaclass adding __contains__ to Enum.

    This allows us to check whether something is a valid alias later in the
    code w/o having to code the try/except everywhere.
    """

    def __contains__(cls, item):
        try:
            cls(item)
        except ValueError:
            return False
        else:
            return True


class RequestType(Enum, metaclass=MemberCheckEnumMeta):
    """A metadata is requested either with a string alias or this enum.

    .. versionadded:: 1.1
    """

    # Metadata is not requested. It will not be routed to the object having the
    # request value as UNREQUESTED.
    UNREQUESTED = False
    # Metadata is requested, and will be routed to the requesting object. There
    # will be no error if the metadata is not provided.
    REQUESTED = True
    # Default metadata request configuration. It should not be passed, and if
    # present, an error is raised for the user to explicitly set the request
    # value.
    ERROR_IF_PASSED = None
    # this is used in `__metadata_request__*` attributes to indicate that a
    # metadata is not present even though it may be present in the
    # corresponding method's signature.
    UNUSED = "$UNUSED$"
    # this is used whenever a default value is changed, and therefore the user
    # should explicitly set the value, otherwise a warning is shown. An example
    # is when a meta-estimator is only a router, but then becomes also a
    # consumer.
    WARN = "$WARN$"


# this is the default used in `set_{method}_request` methods to indicate no change
# requested by the user.
UNCHANGED = "$UNCHANGED$"

# Only the following methods are supported in the routing mechanism. Adding new
# methods at the moment involves monkeypatching this list.
METHODS = [
    "fit",
    "partial_fit",
    "predict",
    "score",
    "split",
    "transform",
    "inverse_transform",
]


# These strings are used to dynamically generate the docstrings for
# set_{method}_request methods.
REQUESTER_DOC = """        Request metadata passed to the ``{method}`` method.

        Please check :ref:`User Guide <metadata_routing>` on how the routing
        mechanism works.

        .. versionadded:: 1.1

        Parameters
        ----------
"""
REQUESTER_DOC_PARAM = """        {metadata} : RequestType, str, True, False, or None, \
                    default=UNCHANGED
            Whether {metadata} should be passed to {method} by meta-estimators or
            not, and if yes, should it have an alias.

            - True or RequestType.REQUESTED: {metadata} is requested, and passed to \
{method} if provided.

            - False or RequestType.UNREQUESTED: {metadata} is not requested and the \
meta-estimator will not pass it to {method}.

            - None or RequestType.ERROR_IF_PASSED: {metadata} is not requested, and \
the meta-estimator will raise an error if the user provides {metadata}

            - str: {metadata} should be passed to the meta-estimator with this given \
alias instead of the original name.

            The default (UNCHANGED) retains the existing request. This allows
            you to change the request for some parameters and not others.

"""
REQUESTER_DOC_RETURN = """        Returns
        -------
        self : estimator
            Returns the object itself.
"""


class MethodMetadataRequest:
    """Contains the metadata request info for a single method.

    Refer to :class:`MetadataRequest` for how this class is used.

    .. versionadded:: 1.1

    Parameters
    ----------
    owner : str
        The name of the object owning these requests.

    method : str
        The name of the method to which these requests belong.
    """

    def __init__(self, owner, method):
        self._requests = dict()
        self.owner = owner
        self.method = method

    @property
    def requests(self):
        return self._requests

    def add_request(
        self,
        *,
        prop,
        alias,
    ):
        """Add request info for a prop.

        Parameters
        ----------
        prop : str
            The property for which a request is set.

        alias : str, RequestType, or {True, False, None}
            The alias which is routed to `prop`

            - str: the name which should be used as an alias when a meta-estimator
              routes the metadata.

            - True or RequestType.REQUESTED: requested

            - False or RequestType.UNREQUESTED: not requested

            - None or RequestType.ERROR_IF_PASSED: error if passed
        """
        try:
            alias = RequestType(alias)
        except ValueError:
            if not (isinstance(alias, str) and alias.isidentifier()):
                raise ValueError(
                    "alias should be either a valid identifier or one of "
                    "{None, True, False}, or a RequestType."
                )

        if alias == prop:
            alias = RequestType.REQUESTED

        if alias == RequestType.UNUSED and prop in self._requests:
            del self._requests[prop]
        else:
            self._requests[prop] = alias

        return self

    def _get_param_names(self, original_names):
        """Get the names of all available metadata.

        This method returns the names of all metadata, even the UNREQUESTED
        ones.

        Parameters
        ----------
        original_names : bool
            Controls whether original or aliased names should be returned. If
            ``True``, aliases are ignored and original names are returned.

        Returns
        -------
        names : set of str
            Returns a set of strings with the names of all parameters.
        """
        return set(
            [
                alias if not original_names and alias not in RequestType else prop
                for prop, alias in self._requests.items()
                if alias not in RequestType
                or RequestType(alias) != RequestType.UNREQUESTED
            ]
        )

    def _check_warnings(self, *, params):
        """Check whether metadata is passed which is marked as WARN.

        If any metadata is passed which is marked as WARN, a warning is raised.

        Parameters
        ----------
        params : dict
            The metadata passed to a method.
        """
        params = {} if params is None else params
        warn_params = {
            prop
            for prop, alias in self._requests.items()
            if alias == RequestType.WARN and prop in params
        }
        for param in warn_params:
            warn(
                f"Support for {param} has recently been added to this class. "
                "To maintain backward compatibility, it is ignored now. "
                "You can set the request value to RequestType.UNREQUESTED "
                "to silence this warning, or to RequestType.REQUESTED to "
                "consume and use the metadata."
            )

    def _route_params(self, params=None):
        """Return the input parameters requested by the method.

        The output of this method can be used directly as the input to the
        corresponding method as extra props.

        Parameters
        ----------
        params : dict
            A dictionary of provided metadata.

        Returns
        -------
        params : Bunch
            A :class:`~utils.Bunch` of {prop: value} which can be given to the
            corresponding method.
        """
        self._check_warnings(params=params)
        params = {} if params is None else params
        args = {arg: value for arg, value in params.items() if value is not None}
        res = Bunch()
        for prop, alias in self._requests.items():
            try:
                alias = RequestType(alias)
            except ValueError:
                # we leave alias as is (str) if it's not a RequestType
                pass

            if alias == RequestType.UNREQUESTED or alias == RequestType.WARN:
                continue
            elif alias == RequestType.REQUESTED and prop in args:
                res[prop] = args[prop]
            elif alias == RequestType.ERROR_IF_PASSED and prop in args:
                raise ValueError(
                    f"{prop} is passed but is not explicitly set as "
                    f"requested or not for {self.owner}.{self.method}"
                )
            elif alias in args:
                res[prop] = args[alias]
        return res

    def _serialize(self):
        """Serialize the object.

        Returns
        -------
        obj : dict
            A serialized version of the instance in the form of a dictionary.
        """
        result = dict()
        # Then parameters with string aliases
        result.update(
            {
                prop: alias
                for prop, alias in self._requests.items()
                if alias not in RequestType
            }
        )
        # And at last the parameters with RequestType routing info
        result.update(
            {
                prop: RequestType(alias)
                for prop, alias in self._requests.items()
                if alias in RequestType
            }
        )
        return result

    def __repr__(self):
        return str(self._serialize())

    def __str__(self):
        return str(repr(self))


class MetadataRequest:
    """Contains the metadata request info of a consumer.

    Instances of :class:`MethodMetadataRequest` are used in this class for each
    available method under `metadatarequest.{method}`.

    Consumers-only classes such as simple estimators return a serialized
    version of this class as the output of `get_metadata_routing()`.

    .. versionadded:: 1.1

    Parameters
    ----------
    owner : str
        The name of the object to which these reqeusts belong.
    """

    # this is here for us to use this attribute's value instead of doing
    # `isinstance`` in our checks, so that we avoid issues when people vendor
    # this file instad of using it directly from scikit-learn.
    type = "request"

    def __init__(self, owner):
        for method in METHODS:
            setattr(self, method, MethodMetadataRequest(owner=owner, method=method))

    def _get_param_names(self, method, original_names, ignore_self=None):
        """Get the names of all available metadata for a method.

        This method returns the names of all metadata, even the UNREQUESTED
        ones.

        Parameters
        ----------
        method : str
            The name of the method for which metadata names are requested.

        original_names : bool
            Controls whether original or aliased names should be returned. If
            ``True``, aliases are ignored and original names are returned.

        ignore_self : bool
            Ignored. Present for API compatibility.

        Returns
        -------
        names : set of str
            Returns a set of strings with the names of all parameters.
        """
        return getattr(self, method)._get_param_names(original_names=original_names)

    def _route_params(self, *, method, params):
        """Return the input parameters requested by a method.

        The output of this method can be used directly as the input to the
        corresponding method as extra props.

        Parameters
        ----------
        method : str
            The name of the method for which the parameters are requested and
            routed.

        params : dict
            A dictionary of provided metadata.

        Returns
        -------
        params : Bunch
            A :class:`~utils.Bunch` of {prop: value} which can be given to the
            corresponding method.
        """
        return getattr(self, method)._route_params(params=params)

    def _check_warnings(self, *, method, params):
        """Check whether metadata is passed which is marked as WARN.

        If any metadata is passed which is marked as WARN, a warning is raised.

        Parameters
        ----------
        method : str
            The name of the method for which the warnings should be checked.

        params : dict
            The metadata passed to a method.
        """
        getattr(self, method)._check_warnings(params=params)

    def _serialize(self):
        """Serialize the object.

        Returns
        -------
        obj : dict
            A serialized version of the instance in the form of a dictionary.
        """
        output = dict()
        for method in METHODS:
            mmr = getattr(self, method)
            if len(mmr.requests):
                output[method] = mmr._serialize()
        return output

    def __repr__(self):
        return str(self._serialize())

    def __str__(self):
        return str(repr(self))


class MethodMapping:
    """Stores the mapping from an object's methods to a router's methods.

    This class is primarily used in a ``get_metadata_routing()`` of a router
    object when defining the mapping between a sub-object (a sub-estimator or a
    scorer) to the router's methods. It stores a collection of ``Route``
    namedtuples.

    .. versionadded:: 1.1
    """

    def __init__(self):
        self._routes = []

    def __iter__(self):
        for route in self._routes:
            yield (route.callee, route.caller)

    def add(self, *, callee, caller):
        """Add a method mapping.

        Parameters
        ----------
        callee : str
            Child object's method name. This method is called in ``caller``.

        caller : str
            Parent estimator's method name in which the ``callee`` is called.

        Returns
        -------
        self : MethodMapping
            Returns self.
        """
        if callee not in METHODS:
            raise ValueError(
                f"Given callee:{callee} is not a valid method. Valid methods are:"
                f" {METHODS}"
            )
        if caller not in METHODS:
            raise ValueError(
                f"Given caller:{caller} is not a valid method. Valid methods are:"
                f" {METHODS}"
            )
        self._routes.append(MethodPair(callee=callee, caller=caller))
        return self

    def _serialize(self):
        """Serialize the object.

        Returns
        -------
        obj : list
            A serialized version of the instance in the form of a list.
        """
        result = list()
        for route in self._routes:
            result.append({"callee": route.callee, "caller": route.caller})
        return result

    @classmethod
    def from_str(cls, route):
        """Construct an instance from a string.

        Parameters
        ----------
        route : str
            A string representing the mapping, it can be:

              - `"one-to-one"`: a one to one mapping for all methods.
              - `"method"`: the name of a single method.

        Returns
        -------
        obj : MethodMapping
            A :class:`~utils.metadata_requests.MethodMapping` instance
            constructed from the given string.
        """
        routing = cls()
        if route == "one-to-one":
            for method in METHODS:
                routing.add(callee=method, caller=method)
        elif route in METHODS:
            routing.add(callee=route, caller=route)
        else:
            raise ValueError("route should be 'one-to-one' or a single method!")
        return routing

    def __repr__(self):
        return str(self._serialize())

    def __str__(self):
        return str(repr(self))


class MetadataRouter:
    """Stores and handles metadata routing for a router object.

    This class is used by router objects to store and handle metadata routing.
    Routing information is stored as a dictionary of the form ``{"object_name":
    RouteMappingPair(method_mapping, routing_info)}``, where ``method_mapping``
    is an instance of :class:`~utils.metadata_requests.MethodMapping` and
    ``routing_info`` is either a
    :class:`~utils.metadata_requests.MetadataRequest` or a
    :class:`~utils.metadata_requests.MetadataRouter` instance.

    .. versionadded:: 1.1

    Parameters
    ----------
    owner : str
        The name of the object to which these requests belong.
    """

    # this is here for us to use this attribute's value instead of doing
    # `isinstance`` in our checks, so that we avoid issues when people vendor
    # this file instad of using it directly from scikit-learn.
    type = "router"

    def __init__(self, owner):
        self._route_mappings = dict()
        self._self = None
        self.owner = owner

    def add_self(self, obj):
        """Add `self` to the routing.

        Parameters
        ----------
        obj : object
            This is typically `self` in a `get_metadata_routing()`.

        Returns
        -------
        self : MetadataRouter
            Returns `self`.
        """
        if getattr(obj, "type", None) == "request":
            self._self = deepcopy(obj)
        elif hasattr(obj, "_get_metadata_request"):
            self._self = deepcopy(obj._get_metadata_request())
        else:
            raise ValueError(
                "Given object is neither a `MetadataRequest` nor does it implement the"
                " require API. Inheriting from `BaseEstimator` implements the required"
                " API."
            )
        return self

    def add(self, *, method_mapping, **objs):
        """Add named objects with their corresponding method mapping.

        Parameters
        ----------
        method_mapping : MethodMapping or str
            The mapping between the child and the parent's methods. If str, the
            output of :func:`~utils.metadata_requests.MethodMapping.from_str`
            is used.

        **objs : dict
            A dictionary of objects from which metadata is extracted by calling
            :func:`~utils.metadata_requests.get_routing_for_object` on them.

        Returns
        -------
        self : MetadataRouter
            Returns `self`.
        """
        if isinstance(method_mapping, str):
            method_mapping = MethodMapping.from_str(method_mapping)
        for name, obj in objs.items():
            self._route_mappings[name] = RouterMappingPair(
                mapping=method_mapping, router=get_routing_for_object(obj)
            )
        return self

    def _get_param_names(self, *, method, original_names, ignore_self):
        """Get the names of all available metadata for a method.

        This method returns the names of all metadata, even the UNREQUESTED
        ones.

        Parameters
        ----------
        method : str
            The name of the method for which metadata names are requested.

        original_names : bool
            Controls whether original or aliased names should be returned,
            which only applies to the stored `self`. If no `self` routing
            object is stored, this parameter has no effect.

        ignore_self : bool
            If `self._self` should be ignored. This is used in
            `_get_squashed_params`. If ``True``, ``original_names`` has no
            effect.

        Returns
        -------
        names : set of str
            Returns a set of strings with the names of all parameters.
        """
        res = set()
        if self._self and not ignore_self:
            res = res.union(
                self._self._get_param_names(
                    method=method, original_names=original_names
                )
            )

        for name, route_mapping in self._route_mappings.items():
            for callee, caller in route_mapping.mapping:
                if caller == method:
                    res = res.union(
                        route_mapping.router._get_param_names(
                            method=callee, original_names=False, ignore_self=False
                        )
                    )
        return set(res)

    def _get_squashed_params(self, *, params, method):
        """Get input for a method of a router w/o validation.

        This is used when a router is used as a child object of another router.
        The parent router then passes all parameters understood by the child
        object to it and delegates their validation to the child.

        The output of this method can be used directly as the input to the
        corresponding method as extra props.

        Parameters
        ----------
        method : str
            The name of the method for which the parameters are requested and
            routed.

        params : dict
            A dictionary of provided metadata.

        Returns
        -------
        params : Bunch
            A :class:`~utils.Bunch` of {prop: value} which can be given to the
            corresponding method.
        """
        res = Bunch()
        if self._self:
            res.update(self._self._route_params(params=params, method=method))

        param_names = self._get_param_names(
            method=method, original_names=False, ignore_self=True
        )
        child_params = {
            key: value for key, value in params.items() if key in param_names
        }
        for key in set(res.keys()).intersection(child_params.keys()):
            # conflicts are okay if the passed objects are the same, but it's
            # an issue if they're different objects.
            if child_params[key] is not res[key]:
                raise ValueError(
                    f"In {self.owner}, there is a conflict on {key} between what is"
                    " requested for this estimator and what is requested by its"
                    " children. You can resolve this conflict by using an alias for"
                    " the child estimator(s) requested metadata."
                )

        res.update(child_params)
        return res

    def route_params(self, *, caller, params):
        """Return the input parameters requested by child objects.

        The output of this method is a bunch, which includes the inputs for all
        methods of each child object that are used in the router's `method`.

        If the router is also a consumer, it also checks for warnings of
        `self`'s/consumer's requested metadata.

        Parameters
        ----------
        caller : str
            The name of the method for which the parameters are requested and
            routed. If called inside the :term:`fit` method of a router, it
            would be `"fit"`.

        params : dict
            A dictionary of provided metadata.

        Returns
        -------
        params : Bunch
            A :class:`~utils.Bunch` of the form
            ``{"object_name": {"method_name": {prop: value}}}`` which can be
            used to pass the required metadata to corresponding methods or
            corresponding child objects.
        """
        if self._self:
            self._self._check_warnings(params=params, method=caller)

        res = Bunch()
        for name, route_mapping in self._route_mappings.items():
            router, mapping = route_mapping.router, route_mapping.mapping

            res[name] = Bunch()
            for _callee, _caller in mapping:
                if _caller == caller:
                    if router.type == "request":
                        res[name][_callee] = router._route_params(
                            params=params, method=_callee
                        )
                    else:
                        res[name][_callee] = router._get_squashed_params(
                            params=params, method=_callee
                        )
        return res

    def validate_metadata(self, *, method, params):
        """Validate given metadata for a method.

        This raises a ``ValueError`` if some of the passed metadata are not
        understood by child objects.

        Parameters
        ----------
        method : str
            The name of the method for which the parameters are requested and
            routed. If called inside the :term:`fit` method of a router, it
            would be `"fit"`.

        params : dict
            A dictionary of provided metadata.
        """
        param_names = self._get_param_names(
            method=method, original_names=True, ignore_self=False
        )
        if self._self:
            self_params = self._self._get_param_names(
                method=method, original_names=True
            )
        else:
            self_params = set()
        extra_keys = set(params.keys()) - param_names - self_params
        if extra_keys:
            raise TypeError(
                f"{method} got unexpected argument(s) {extra_keys}, which are "
                "not requested metadata in any object."
            )

    def _serialize(self):
        """Serialize the object.

        Returns
        -------
        obj : dict
            A serialized version of the instance in the form of a dictionary.
        """
        res = dict()
        if self._self:
            res["$self"] = self._self._serialize()
        for name, route_mapping in self._route_mappings.items():
            res[name] = dict()
            res[name]["mapping"] = route_mapping.mapping._serialize()
            res[name]["router"] = route_mapping.router._serialize()

        return res

    def __iter__(self):
        if self._self:
            yield "$self", RouterMappingPair(
                mapping=MethodMapping.from_str("one-to-one"), router=self._self
            )
        for name, route_mapping in self._route_mappings.items():
            yield (name, route_mapping)

    def __repr__(self):
        return str(self._serialize())

    def __str__(self):
        return str(repr(self))


def get_routing_for_object(obj=None):
    """Get a ``Metadata{Router, Request}`` instance from the given object.

    This factory function returns a
    :class:`~utils.metadata_request.MetadataRouter` or a
    :class:`~utils.metadata_request.MetadataRequest` from the given input.

    This function always returns a copy or an instance constructed from the
    intput, such that changing the output of this function will not change the
    original object.

    .. versionadded:: 1.1

    Parameters
    ----------
    obj : object
        - If the object is already a
            :class:`~utils.metadata_requests.MetadataRequest` or a
            :class:`~utils.metadata_requests.MetadataRouter`, return a copy
            of that.
        - If the object provides a `get_metadata_routing` method, return a copy
            of the output of the method.
        - If the object is a dict, return the deserialized instance from it.
        - If ``None``, return an empty
            :class:`~utils.metadata_requests.MetadataRequest`.

    Returns
    -------
    obj : MetadataRequest or MetadataRouting
        A ``MetadataRequest`` or a ``MetadataRouting`` taken or created from
        the given object.
    """
    if obj is None:
        return MetadataRequest(owner=None)

    # doing this instead of a try/except since an AttributeError could be raised
    # for other reasons.
    if hasattr(obj, "get_metadata_routing"):
        return deepcopy(obj.get_metadata_routing())

    if getattr(obj, "type", None) in ["request", "router"]:
        return deepcopy(obj)

    return MetadataRequest(owner=None)


class RequestMethod:
    """
    A descriptor for request methods.

    .. versionadded:: 1.1

    Parameters
    ----------
    name : str
        The name of the method for which the request function should be
        created, e.g. ``"fit"`` would create a ``set_fit_request`` function.

    keys : list of str
        A list of strings which are accepted parameters by the created
        function, e.g. ``["sample_weight"]`` if the corresponding method
        accepts it as a metadata.

    Notes
    -----
    This class is a descriptor [1]_ and uses PEP-362 to set the signature of
    the returned function [2]_.

    References
    ----------
    .. [1] https://docs.python.org/3/howto/descriptor.html

    .. [2] https://www.python.org/dev/peps/pep-0362/
    """

    def __init__(self, name, keys):
        self.name = name
        self.keys = keys

    def __get__(self, instance, owner):
        # we would want to have a method which accepts only the expected args
        def func(**kw):
            if set(kw) - set(self.keys):
                raise TypeError(
                    f"Unexpected args: {set(kw) - set(self.keys)}. Accepted arguments"
                    f" are: {set(self.keys)}"
                )

            requests = instance._get_metadata_request()
            method_metadata_request = getattr(requests, self.name)

            for prop, alias in kw.items():
                if alias is not UNCHANGED:
                    method_metadata_request.add_request(prop=prop, alias=alias)
            instance._metadata_request = requests

            return instance

        # Now we set the relevant attributes of the function so that it seems
        # like a normal method to the end user, with known expected arguments.
        func.__name__ = f"set_{self.name}_request"
        params = [
            inspect.Parameter(
                name="self",
                kind=inspect.Parameter.POSITIONAL_OR_KEYWORD,
                annotation=owner,
            )
        ]
        params.extend(
            [
                inspect.Parameter(
                    k,
                    inspect.Parameter.KEYWORD_ONLY,
                    default=UNCHANGED,
                    annotation=Optional[Union[RequestType, str]],
                )
                for k in self.keys
            ]
        )
        func.__signature__ = inspect.Signature(
            params,
            return_annotation=owner,
        )
        doc = REQUESTER_DOC.format(method=self.name)
        for metadata in self.keys:
            doc += REQUESTER_DOC_PARAM.format(metadata=metadata, method=self.name)
        doc += REQUESTER_DOC_RETURN
        func.__doc__ = doc
        return func


class _MetadataRequester:
    """Mixin class for adding metadata request functionality.

    .. versionadded:: 1.1
    """

    def __init_subclass__(cls, **kwargs):
        """Set the ``set_{method}_request`` methods.

        This uses PEP-487 [1]_ to set the ``set_{method}_request`` methods. It
        looks for the information available in the set default values which are
        set using ``__metadata_request__*`` class attributes.

        References
        ----------
        .. [1] https://www.python.org/dev/peps/pep-0487
        """
        try:
            requests = cls._get_default_requests()
        except Exception:
            # if there are any issues in the default values, it will be raised
            # when ``get_metadata_routing`` is called. Here we are going to
            # ignore all the issues such as bad defaults etc.
            super().__init_subclass__(**kwargs)
            return

        for method in METHODS:
            mmr = getattr(requests, method)
            # set ``set_{method}_request``` methods
            if not len(mmr.requests):
                continue
            setattr(
                cls,
                f"set_{method}_request",
                RequestMethod(method, sorted(mmr.requests.keys())),
            )
        super().__init_subclass__(**kwargs)

    @classmethod
    def _get_default_requests(cls):
        """Collect default request values.

        This method combines the information present in ``metadata_request__*``
        class attributes.
        """

        requests = MetadataRequest(owner=cls.__name__)

        # need to go through the MRO since this is a class attribute and
        # ``vars`` doesn't report the parent class attributes. We go through
        # the reverse of the MRO since cls is the first in the tuple and object
        # is the last.
        defaults = dict()
        for klass in reversed(inspect.getmro(cls)):
            klass_defaults = {
                attr: value
                for attr, value in vars(klass).items()
                if "__metadata_request__" in attr
            }
            defaults.update(klass_defaults)
        defaults = dict(sorted(defaults.items()))

        # First take all arguments from the method signatures and have them as
        # ERROR_IF_PASSED, except X, y, *args, and **kwargs.
        for method in METHODS:
            # Here we use `isfunction` instead of `ismethod` because calling `getattr`
            # on a class instead of an instance returns an unbound function.
            if not hasattr(cls, method) or not inspect.isfunction(getattr(cls, method)):
                continue
            # ignore the first parameter of the method, which is usually "self"
            params = list(inspect.signature(getattr(cls, method)).parameters.items())[
                1:
            ]
            for pname, param in params:
                if pname in {"X", "y", "Y"}:
                    continue
                if param.kind in {param.VAR_POSITIONAL, param.VAR_KEYWORD}:
                    continue
                getattr(requests, method).add_request(
                    prop=pname,
                    alias=RequestType.ERROR_IF_PASSED,
                )

        # Then overwrite those defaults with the ones provided in
        # __metadata_request__* attributes, which are provided in `requests` here.

        for attr, value in defaults.items():
            # we don't check for attr.startswith() since python prefixes attrs
            # starting with __ with the `_ClassName`.
            substr = "__metadata_request__"
            method = attr[attr.index(substr) + len(substr) :]
            for prop, alias in value.items():
                getattr(requests, method).add_request(prop=prop, alias=alias)
        return requests

    def _get_metadata_request(self):
        """Get requested data properties.

        Please check :ref:`User Guide <metadata_routing>` on how the routing
        mechanism works.

        Returns
        -------
        request : MetadataRequest
            A :class:`~.utils.metadata_requests.MetadataRequest` instance.
        """
        if hasattr(self, "_metadata_request"):
            requests = get_routing_for_object(self._metadata_request)
        else:
            requests = self._get_default_requests()

        return requests

    def get_metadata_routing(self):
        """Get metadata routing of this object.

        Please check :ref:`User Guide <metadata_routing>` on how the routing
        mechanism works.

        Returns
        -------
        routing : dict
            A dict representing a MetadataRouter.
        """
        return self._get_metadata_request()


def process_routing(obj, method, other_params, **kwargs):
    """Validate and route input parameters.

    This function is used inside a router's method, e.g. :term:`fit`,
    to validate the metadata and handle the routing.

    Assuming this signature: ``fit(self, X, y, sample_weight=None, **fit_params)``,
    a call to this function would be:
    ``process_routing(self, fit_params, sample_weight=sample_weight)``.

    .. versionadded:: 1.1

    Parameters
    ----------
    obj : object
        An object implementing ``get_metadata_routing``. Typically a
        meta-estimator.

    method : str
        The name of the router's method in which this function is called.

    other_params : dict
        A dictionary of extra parameters passed to the router's method,
        e.g. ``**fit_params`` passed to a meta-estimator's :term:`fit`.

    **kwargs : dict
        Parameters explicitly accepted and included in the router's method
        signature.

    Returns
    -------
    routed_params : Bunch
        A :class:`~utils.Bunch` with a key for each object added in
        ``get_metadata_routing``, and a key for each method of the child object
        which is used in the router's `method`.
    """
    if not hasattr(obj, "get_metadata_routing"):
        raise AttributeError(
            f"This {repr(obj.__class__.__name__)} has not implemented the routing"
            " method `get_metadata_routing`."
        )
    if method not in METHODS:
        raise TypeError(
            f"Can only route and process input on these methods: {METHODS}, "
            f"while the passed method is: {method}."
        )

    # We take the extra params (**fit_params) which is passed as `other_params`
    # and add the explicitly passed parameters (passed as **kwargs) to it. This
    # is equivalent to a code such as this in a router:
    # if sample_weight is not None:
    #     fit_params["sample_weight"] = sample_weight
    all_params = other_params if other_params is not None else dict()
    all_params.update(kwargs)

    request_routing = get_routing_for_object(obj)
    request_routing.validate_metadata(params=all_params, method=method)
    routed_params = request_routing.route_params(params=all_params, caller=method)

    return routed_params
