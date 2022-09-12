"""
Metadata Routing Utility
"""

# Author: Adrin Jalali <adrin.jalali@gmail.com>
# License: BSD 3 clause

import inspect
from copy import deepcopy
from enum import Enum
from warnings import warn
from collections import namedtuple
from typing import Union, Optional
from ._bunch import Bunch
from ..exceptions import UnsetMetadataPassedError

# This namedtuple is used to store a (mapping, routing) pair. Mapping is a
# MethodMapping object, and routing is the output of `get_metadata_routing`.
# MetadataRouter stores a collection of these namedtuples.
RouterMappingPair = namedtuple("RouterMappingPair", ["mapping", "router"])

# A namedtuple storing a single method route. A collection of these namedtuples
# is stored in a MetadataRouter.
MethodPair = namedtuple("MethodPair", ["callee", "caller"])


class RequestType(Enum):
    """A metadata is requested either with a string alias or this enum.

    .. versionadded:: 1.2
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
    # consumer in a new release.
    WARN = "$WARN$"

    @classmethod
    def is_alias(cls, item):
        """Check if an item is a valid alias.

        Parameters
        ----------
        item : object
            The given item to be checked if it can be an alias.

        Returns
        -------
        result : bool
            Whether the given item is a valid alias.
        """
        try:
            cls(item)
        except ValueError:
            # item is only an alias if it's a valid identifier
            return isinstance(item, str) and item.isidentifier()
        else:
            return False

    @classmethod
    def is_valid(cls, item):
        """Check if an item is a valid RequestType (and not an alias).

        Parameters
        ----------
        item : object
            The given item to be checked.

        Returns
        -------
        result : bool
            Whether the given item is valid.
        """
        try:
            cls(item)
            return True
        except ValueError:
            return False


# this is the default used in `set_{method}_request` methods to indicate no change
# requested by the user.
UNCHANGED = "$UNCHANGED$"

# Only the following methods are supported in the routing mechanism. Adding new
# methods at the moment involves monkeypatching this list.
METHODS = [
    "fit",
    "partial_fit",
    "predict",
    "predict_proba",
    "predict_log_proba",
    "decision_function",
    "score",
    "split",
    "transform",
    "inverse_transform",
]


# These strings are used to dynamically generate the docstrings for
# set_{method}_request methods.
REQUESTER_DOC = """        Request metadata passed to the ``{method}`` method.

        Please see :ref:`User Guide <metadata_routing>` on how the routing
        mechanism works.

        .. versionadded:: 1.2

        Parameters
        ----------
"""
REQUESTER_DOC_PARAM = """        {metadata} : RequestType, str, True, False, or None, \
                    default=UNCHANGED
            Whether ``{metadata}`` should be passed to ``{method}`` by
            meta-estimators or not, and if yes, should it have an alias.

            - True or RequestType.REQUESTED: ``{metadata}`` is requested, and \
passed to ``{method}`` if provided. The request is ignored if \
``{metadata}`` is not provided.

            - False or RequestType.UNREQUESTED: ``{metadata}`` is not requested \
and the meta-estimator will not pass it to ``{method}``.

            - None or RequestType.ERROR_IF_PASSED: ``{metadata}`` is not \
requested, and the meta-estimator will raise an error if the user provides \
``{metadata}``.

            - str: ``{metadata}`` should be passed to the meta-estimator with \
this given alias instead of the original name.

            The default (UNCHANGED) retains the existing request. This allows
            you to change the request for some parameters and not others.

"""
REQUESTER_DOC_RETURN = """        Returns
        -------
        self : object
            The updated object.
"""


class MethodMetadataRequest:
    """A prescription of how metadata is to be passed to a single method.

    Refer to :class:`MetadataRequest` for how this class is used.

    .. versionadded:: 1.2

    Parameters
    ----------
    owner : str
        A display name for the object owning these requests.

    method : str
        The name of the method to which these requests belong.
    """

    def __init__(self, router, owner, method):
        self._requests = dict()
        self.router = router
        self.owner = owner
        self.method = method

    @property
    def requests(self):
        """Dictionary of the form: ``{key: alias}``."""
        return self._requests

    def add_request(
        self,
        *,
        param,
        alias,
    ):
        """Add request info for a metadata.

        Parameters
        ----------
        param : str
            The property for which a request is set.

        alias : str, RequestType, or {True, False, None}
            Specifies which metadata should be routed to `param`

            - str: the name (or alias) of metadata given to a meta-estimator that
              should be routed to this parameter.

            - True or RequestType.REQUESTED: requested

            - False or RequestType.UNREQUESTED: not requested

            - None or RequestType.ERROR_IF_PASSED: error if passed
        """
        if RequestType.is_valid(alias):
            alias = RequestType(alias)
        elif not RequestType.is_alias(alias):
            raise ValueError(
                "alias should be either a valid identifier or one of "
                "{None, True, False}, or a RequestType."
            )

        if alias != self._requests.get(param, None):
            self.router._is_default_request = False

        if alias == param:
            alias = RequestType.REQUESTED

        if alias == RequestType.UNUSED and param in self._requests:
            del self._requests[param]
        else:
            self._requests[param] = alias

        return self

    def _get_param_names(self, return_alias):
        """Get names of all metadata that can be consumed or routed by this method.

        This method returns the names of all metadata, even the UNREQUESTED
        ones.

        Parameters
        ----------
        return_alias : bool
            Controls whether original or aliased names should be returned. If
            ``False``, aliases are ignored and original names are returned.

        Returns
        -------
        names : set of str
            A set of strings with the names of all parameters.
        """
        return set(
            alias if return_alias and not RequestType.is_valid(alias) else prop
            for prop, alias in self._requests.items()
            if not RequestType.is_valid(alias)
            or RequestType(alias) != RequestType.UNREQUESTED
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
        """Prepare the given parameters to be passed to the method.

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
        unrequested = dict()
        params = {} if params is None else params
        args = {arg: value for arg, value in params.items() if value is not None}
        res = Bunch()
        for prop, alias in self._requests.items():
            if RequestType.is_valid(alias):
                alias = RequestType(alias)

            if alias == RequestType.UNREQUESTED or alias == RequestType.WARN:
                continue
            elif alias == RequestType.REQUESTED and prop in args:
                res[prop] = args[prop]
            elif alias == RequestType.ERROR_IF_PASSED and prop in args:
                unrequested[prop] = args[prop]
            elif alias in args:
                res[prop] = args[alias]
        if unrequested:
            raise UnsetMetadataPassedError(
                message=(
                    f"[{', '.join([key for key in unrequested])}] are passed but are"
                    " not explicitly set as requested or not for"
                    f" {self.owner}.{self.method}"
                ),
                unrequested_params=unrequested,
                routed_params=res,
            )
        return res

    def _serialize(self):
        """Serialize the object.

        Returns
        -------
        obj : dict
            A serialized version of the instance in the form of a dictionary.
        """
        return {
            prop: RequestType(alias) if RequestType.is_valid(alias) else alias
            for prop, alias in self._requests.items()
        }

    def __repr__(self):
        return str(self._serialize())

    def __str__(self):
        return str(repr(self))


class MetadataRequest:
    """Contains the metadata request info of a consumer.

    Instances of :class:`MethodMetadataRequest` are used in this class for each
    available method under `metadatarequest.{method}`.

    Consumer-only classes such as simple estimators return a serialized
    version of this class as the output of `get_metadata_routing()`.

    .. versionadded:: 1.2

    Parameters
    ----------
    owner : str
        The name of the object to which these requests belong.
    """

    # this is here for us to use this attribute's value instead of doing
    # `isinstance` in our checks, so that we avoid issues when people vendor
    # this file instead of using it directly from scikit-learn.
    _type = "metadata_request"

    def __init__(self, owner):
        # this is used to check if the user has set any request values
        self._is_default_request = False
        for method in METHODS:
            setattr(
                self,
                method,
                MethodMetadataRequest(router=self, owner=owner, method=method),
            )

    def _get_param_names(self, method, return_alias, ignore_self=None):
        """Get names of all metadata that can be consumed or routed by specified \
            method.

        This method returns the names of all metadata, even the UNREQUESTED
        ones.

        Parameters
        ----------
        method : str
            The name of the method for which metadata names are requested.

        return_alias : bool
            Controls whether original or aliased names should be returned. If
            ``False``, aliases are ignored and original names are returned.

        ignore_self : bool
            Ignored. Present for API compatibility.

        Returns
        -------
        names : set of str
            A set of strings with the names of all parameters.
        """
        return getattr(self, method)._get_param_names(return_alias=return_alias)

    def _route_params(self, *, method, params):
        """Prepare the given parameters to be passed to the method.

        The output of this method can be used directly as the input to the
        corresponding method as extra keyword arguments to pass metadata.

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
    """Stores the mapping between callee and caller methods for a router.

    This class is primarily used in a ``get_metadata_routing()`` of a router
    object when defining the mapping between a sub-object (a sub-estimator or a
    scorer) to the router's methods. It stores a collection of ``Route``
    namedtuples.

    Iterating through an instance of this class will yield named
    ``MethodPair(callee, caller)`` tuples.

    .. versionadded:: 1.2
    """

    def __init__(self):
        self._routes = []

    def __iter__(self):
        return iter(self._routes)

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

    .. versionadded:: 1.2

    Parameters
    ----------
    owner : str
        The name of the object to which these requests belong.
    """

    # this is here for us to use this attribute's value instead of doing
    # `isinstance`` in our checks, so that we avoid issues when people vendor
    # this file instad of using it directly from scikit-learn.
    _type = "metadata_router"

    def __init__(self, owner):
        self._route_mappings = dict()
        # `_self` is used if the router is also a consumer. _self, (added using
        # `add_self()`) is treated differently from the other objects which are
        # stored in _route_mappings.
        self._self = None
        # this attribute is used to decide if there should be an error raised
        # or a FutureWarning if a metadata is passed which is not requested.
        self._warn_on = dict()
        self.owner = owner

    @property
    def _is_default_request(self):
        """Return ``True`` only if all sub-components have default values."""
        if self._self and not self._self._is_default_request:
            return False

        for router_mapping in self._route_mappings.values():
            if not router_mapping.router._is_default_request:
                return False

        return True

    def add_self(self, obj):
        """Add `self` (as a consumer) to the routing.

        This method is used if the router is also a consumer, and hence the
        router itself needs to be included in the routing. The passed object
        can be an estimator or a
        :class:``~utils.metadata_requests.MetadataRequest``.

        A router should add itself using this method instead of `add` since it
        should be treated differently than the other objects to which metadata
        is routed by the router.

        Parameters
        ----------
        obj : object
            This is typically the router instance, i.e. `self` in a
            ``get_metadata_routing()`` implementation. It can also be a
            ``MetadataRequest`` instance.

        Returns
        -------
        self : MetadataRouter
            Returns `self`.
        """
        if getattr(obj, "_type", None) == "metadata_request":
            self._self = deepcopy(obj)
        elif hasattr(obj, "_get_metadata_request"):
            self._self = deepcopy(obj._get_metadata_request())
        else:
            raise ValueError(
                "Given `obj` is neither a `MetadataRequest` nor does it implement the"
                " required API. Inheriting from `BaseEstimator` implements the required"
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
        else:
            method_mapping = deepcopy(method_mapping)

        for name, obj in objs.items():
            self._route_mappings[name] = RouterMappingPair(
                mapping=method_mapping, router=get_routing_for_object(obj)
            )
        return self

    def _get_param_names(self, *, method, return_alias, ignore_self):
        """Get names of all metadata that can be consumed or routed by specified \
            method.

        This method returns the names of all metadata, even the UNREQUESTED
        ones.

        Parameters
        ----------
        method : str
            The name of the method for which metadata names are requested.

        return_alias : bool
            Controls whether original or aliased names should be returned,
            which only applies to the stored `self`. If no `self` routing
            object is stored, this parameter has no effect.

        ignore_self : bool
            If `self._self` should be ignored. This is used in `_route_params`.
            If ``True``, ``return_alias`` has no effect.

        Returns
        -------
        names : set of str
            A set of strings with the names of all parameters.
        """
        res = set()
        if self._self and not ignore_self:
            res = res.union(
                self._self._get_param_names(method=method, return_alias=return_alias)
            )

        for name, route_mapping in self._route_mappings.items():
            for callee, caller in route_mapping.mapping:
                if caller == method:
                    res = res.union(
                        route_mapping.router._get_param_names(
                            method=callee, return_alias=True, ignore_self=False
                        )
                    )
        return set(res)

    def _route_params(self, *, params, method):
        """Prepare the given parameters to be passed to the method.

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
            method=method, return_alias=True, ignore_self=True
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
        methods of each child object that are used in the router's `caller`
        method.

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
                    res[name][_callee] = self._route_warn_or_error(
                        child=name, router=router, params=params, method=_callee
                    )
        return res

    def _route_warn_or_error(self, *, child, router, params, method):
        """Route parameters while handling error or deprecation warning choice.

        This method warns instead of raising an error if the parent object
        has set ``warn_on`` for the child object's method and the user has not
        set any metadata request for that child object. This is used during the
        deprecation cycle for backward compatibility.

        Parameters
        ----------
        child : str
            The name of the child object.

        router : MetadataRouter or MetadataRequest
            The router for the child object.

        params : dict
            The parameters to be routed.

        method : str
            The name of the callee method.

        Returns
        -------
        dict
            The routed parameters.
        """
        try:
            routed_params = router._route_params(params=params, method=method)
        except UnsetMetadataPassedError as e:
            warn_on = self._warn_on.get(child, {})
            if method not in warn_on:
                # there is no warn_on set for this method of this child object,
                # we raise as usual.
                raise
            if not router._is_default_request:
                # the user has set at least one request value for this child
                # object, but not for all of them. Therefore we raise as usual.
                raise
            # now we move everything which has a warn_on flag from
            # `unrequested_params` to routed_params, and then raise if anything
            # is left. Otherwise we have a perfectly formed `routed_params` and
            # we return that.
            warn_on_params = warn_on.get(method, {"params": [], "raise_on": "1.4"})
            warn_keys = list(e.unrequested_params.keys())
            routed_params = e.routed_params
            # if params is None, we accept and warn on everything.
            warn_params = warn_on_params["params"]
            if warn_params is None:
                warn_params = warn_keys

            for param in warn_params:
                if param in e.unrequested_params:
                    routed_params[param] = e.unrequested_params.pop(param)

            # check if anything is left, and if yes, we raise as usual
            if e.unrequested_params:
                raise

            # Finally warn before returning the routed parameters.
            warn(
                "You are passing metadata for which the request values are not"
                f" explicitly set: {', '.join(warn_keys)}. From version"
                f" {warn_on_params['raise_on']} this results in the following error:"
                f" {str(e)}",
                FutureWarning,
            )
        return routed_params

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
            method=method, return_alias=False, ignore_self=False
        )
        if self._self:
            self_params = self._self._get_param_names(method=method, return_alias=False)
        else:
            self_params = set()
        extra_keys = set(params.keys()) - param_names - self_params
        if extra_keys:
            raise TypeError(
                f"{method} got unexpected argument(s) {extra_keys}, which are "
                "not requested metadata in any object."
            )

    def warn_on(self, *, child, method, params, raise_on="1.4"):
        """Set where deprecation warnings on no set requests should occur.

        This method is used in meta-estimators during the transition period for
        backward compatibility. Expected behavior for meta-estimators on a code
        such as ``RFE(Ridge()).fit(X, y, sample_weight=sample_weight)`` is to
        raise a ``ValueError`` complaining about the fact that ``Ridge()`` has
        not explicitly set the request values for ``sample_weight``. However,
        this breaks backward compatibility for existing meta-estimators.

        Calling this method on a ``MetadataRouter`` object such as
        ``warn_on(child='estimator', method='fit', params=['sample_weight'])``
        tells the router to raise a ``FutureWarning`` instead of a
        ``ValueError`` if the child object has no set requests for
        ``sample_weight`` during ``fit``.

        You can find more information on how to use this method in the
        developer guide:
        :ref:`sphx_glr_auto_examples_plot_metadata_routing.py`.

        Parameters
        ----------
        child : str
            The name of the child object. The names come from the keyword
            arguments passed to the ``add`` method.

        method : str
            The method for which there should be a ``FutureWarning``
            instead of a ``ValueError`` for given params.

        params : list of str
            The list of parameters for which there should be a
            ``FutureWarning`` instead of a ``ValueError``. If ``None``, the
            rule is applied on all parameters.

        raise_on : str, default="1.4"
            The version after which there should be an error. Used in the
            warning message to inform users.

        Returns
        -------
        self : MetadataRouter
            Returns `self`.
        """
        if child not in self._route_mappings:
            raise ValueError(f"Unknown child object: {child}")

        if child not in self._warn_on:
            self._warn_on[child] = dict()
        self._warn_on[child][method] = {"params": params, "raise_on": raise_on}
        return self

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

    This function returns a
    :class:`~utils.metadata_request.MetadataRouter` or a
    :class:`~utils.metadata_request.MetadataRequest` from the given input.

    This function always returns a copy or an instance constructed from the
    intput, such that changing the output of this function will not change the
    original object.

    .. versionadded:: 1.2

    Parameters
    ----------
    obj : object
        - If the object is already a
            :class:`~utils.metadata_requests.MetadataRequest` or a
            :class:`~utils.metadata_requests.MetadataRouter`, return a copy
            of that.
        - If the object provides a `get_metadata_routing` method, return a copy
            of the output of that method.
        - Returns an empty :class:`~utils.metadata_requests.MetadataRequest`
            otherwise.

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

    if getattr(obj, "_type", None) in ["metadata_request", "metadata_router"]:
        return deepcopy(obj)

    return MetadataRequest(owner=None)


class RequestMethod:
    """
    A descriptor for request methods.

    .. versionadded:: 1.2

    Parameters
    ----------
    name : str
        The name of the method for which the request function should be
        created, e.g. ``"fit"`` would create a ``set_fit_request`` function.

    keys : list of str
        A list of strings which are accepted parameters by the created
        function, e.g. ``["sample_weight"]`` if the corresponding method
        accepts it as a metadata.

    validate_keys : bool, default=True
        Whether to check if the requested parameters fit the actual parameters
        of the method.

    Notes
    -----
    This class is a descriptor [1]_ and uses PEP-362 to set the signature of
    the returned function [2]_.

    References
    ----------
    .. [1] https://docs.python.org/3/howto/descriptor.html

    .. [2] https://www.python.org/dev/peps/pep-0362/
    """

    def __init__(self, name, keys, validate_keys=True):
        self.name = name
        self.keys = keys
        self.validate_keys = validate_keys

    def __get__(self, instance, owner):
        # we would want to have a method which accepts only the expected args
        def func(**kw):
            """Updates the request for provided parameters

            This docstring is overwritten below.
            See REQUESTER_DOC for expected functionality
            """
            if self.validate_keys and (set(kw) - set(self.keys)):
                raise TypeError(
                    f"Unexpected args: {set(kw) - set(self.keys)}. Accepted arguments"
                    f" are: {set(self.keys)}"
                )

            requests = instance._get_metadata_request()
            method_metadata_request = getattr(requests, self.name)

            for prop, alias in kw.items():
                if alias is not UNCHANGED:
                    method_metadata_request.add_request(param=prop, alias=alias)
                    requests._is_default_request = False
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

    .. versionadded:: 1.2
    """

    def __init_subclass__(cls, **kwargs):
        """Set the ``set_{method}_request`` methods.

        This uses PEP-487 [1]_ to set the ``set_{method}_request`` methods. It
        looks for the information available in the set default values which are
        set using ``__metadata_request__*`` class attributes, or inferred
        from method signatures.

        The ``__metadata_request__*`` class attributes are used when a method
        does not explicitly accept a metadata through its arguments or if the
        developer would like to specify a request value for those metadata
        which are different from the default ``RequestType.ERROR_IF_PASSED``.

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
    def _build_request_for_signature(cls, router, method):
        """Build the `MethodMetadataRequest` for a method using its signature.

        This method takes all arguments from the method signature and uses
        ``RequestType.ERROR_IF_PASSED`` as their default request value, except
        ``X``, ``y``, ``*args``, and ``**kwargs``.

        Parameters
        ----------
        router : MetadataRequest
            The parent object for the created `MethodMetadataRequest`.
        method : str
            The name of the method.

        Returns
        -------
        method_request : MethodMetadataRequest
            The prepared request using the method's signature.
        """
        mmr = MethodMetadataRequest(router=router, owner=cls.__name__, method=method)
        # Here we use `isfunction` instead of `ismethod` because calling `getattr`
        # on a class instead of an instance returns an unbound function.
        if not hasattr(cls, method) or not inspect.isfunction(getattr(cls, method)):
            return mmr
        # ignore the first parameter of the method, which is usually "self"
        params = list(inspect.signature(getattr(cls, method)).parameters.items())[1:]
        for pname, param in params:
            if pname in {"X", "y", "Y"}:
                continue
            if param.kind in {param.VAR_POSITIONAL, param.VAR_KEYWORD}:
                continue
            mmr.add_request(
                param=pname,
                alias=RequestType.ERROR_IF_PASSED,
            )
        return mmr

    @classmethod
    def _get_default_requests(cls):
        """Collect default request values.

        This method combines the information present in ``metadata_request__*``
        class attributes, as well as determining request keys from method
        signatures.
        """
        requests = MetadataRequest(owner=cls.__name__)

        for method in METHODS:
            setattr(
                requests,
                method,
                cls._build_request_for_signature(router=requests, method=method),
            )

        # Then overwrite those defaults with the ones provided in
        # __metadata_request__* attributes. Defaults set in
        # __metadata_request__* attributes take precedence over signature
        # sniffing.

        # need to go through the MRO since this is a class attribute and
        # ``vars`` doesn't report the parent class attributes. We go through
        # the reverse of the MRO so that child classes have precedence over
        # their parents.
        defaults = dict()
        for base_class in reversed(inspect.getmro(cls)):
            base_defaults = {
                attr: value
                for attr, value in vars(base_class).items()
                if "__metadata_request__" in attr
            }
            defaults.update(base_defaults)
        defaults = dict(sorted(defaults.items()))

        for attr, value in defaults.items():
            # we don't check for attr.startswith() since python prefixes attrs
            # starting with __ with the `_ClassName`.
            substr = "__metadata_request__"
            method = attr[attr.index(substr) + len(substr) :]
            for prop, alias in value.items():
                getattr(requests, method).add_request(param=prop, alias=alias)

        # this indicates that the user has not set any request values for this
        # object
        requests._is_default_request = True

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
        routing : MetadataRequest
            A :class:`~utils.metadata_routing.MetadataRequest` encapsulating
            routing information.
        """
        return self._get_metadata_request()


def process_routing(obj, method, other_params, **kwargs):
    """Validate and route input parameters.

    This function is used inside a router's method, e.g. :term:`fit`,
    to validate the metadata and handle the routing.

    Assuming this signature: ``fit(self, X, y, sample_weight=None, **fit_params)``,
    a call to this function would be:
    ``process_routing(self, fit_params, sample_weight=sample_weight)``.

    .. versionadded:: 1.2

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
        A :class:`~utils.Bunch` of the form ``{"object_name": {"method_name":
        {prop: value}}}`` which can be used to pass the required metadata to
        corresponding methods or corresponding child objects. The object names
        are those defined in `obj.get_metadata_routing()`.
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
    all_params = {
        param: value for param, value in all_params.items() if value is not None
    }

    request_routing = get_routing_for_object(obj)
    request_routing.validate_metadata(params=all_params, method=method)
    routed_params = request_routing.route_params(params=all_params, caller=method)

    return routed_params
