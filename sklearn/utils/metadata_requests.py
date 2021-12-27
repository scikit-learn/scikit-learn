"""
Metadata Routing Utility
"""

# Author: Adrin Jalali <adrin.jalali@gmail.com>
# License: BSD 3 clause

import functools
import inspect
from copy import deepcopy
from enum import Enum
from warnings import warn
from collections import namedtuple
from typing import Union, Optional
from ..externals._sentinels import sentinel  # type: ignore # mypy error!!!
from ._bunch import Bunch

# This namedtuple is used to store a (mapping, routing) pair. Mapping is a
# MethodMapping object, and routing is the output of `get_metadata_routing`.
# MetadataRouter stores a collection of these namedtuples.
RouteMappingPair = namedtuple("RouteMappingPair", ["mapping", "routing"])

# A namedtuple storing a single method route. A collection of these namedtuples
# is stored in a MetadataRouter.
Route = namedtuple("Route", ["method", "used_in"])


class RequestType(Enum):
    """A metadata is requested either with a string alias or this enum.

    .. versionadded: 1.1
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
    # this sentinel is used in `__metadata_request__*` attributes to indicate
    # that a metadata is not present even though it may be present in the
    # corresponding method's signature.
    UNUSED = sentinel("UNUSED")
    # this sentinel is used whenever a default value is changed, and therefore
    # the user should explicitly set the value, otherwise a warning is shown.
    # An example is when a meta-estimator is only a router, but then becomes
    # also a consumer.
    WARN = sentinel("WARN")


# this sentinel is the default used in `{method}_requests` methods to indicate
# no change requested by the user.
UNCHANGED = sentinel("UNCHANGED")

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


REQUESTER_DOC = """        Request metadata passed to the ``{method}`` method.

        Please check :ref:`User Guide <metadata_routing>` on how the routing
        mechanism works.

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

"""
REQUESTER_DOC_RETURN = """        Returns
        -------
        self
            Returns the object itself.
"""


class MethodMetadataRequest:
    """Contains the metadata request info for a single method.

    Refer to :class:`MetadataRequest` for how this class is used.

    .. versionadded:: 1.1

    Parameters
    ----------
    name : str
        The name of the method to which these requests belong.
    """

    def __init__(self, name):
        self._requests = dict()
        self.name = name

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
        if not isinstance(alias, str):
            try:
                alias = RequestType(alias)
            except ValueError:
                raise ValueError(
                    "alias should be either a string or one of "
                    "{None, True, False}, or a RequestType."
                )

        if alias == prop:
            alias = RequestType.REQUESTED

        if alias == RequestType.UNUSED and prop in self._requests:
            del self._requests[prop]
        else:
            self._requests[prop] = alias

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
            sorted(
                [
                    alias if not original_names and isinstance(alias, str) else prop
                    for prop, alias in self._requests.items()
                    if isinstance(alias, str)
                    or RequestType(alias) != RequestType.UNREQUESTED
                ]
            )
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

    def _get_params(self, params=None):
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
            if not isinstance(alias, str):
                alias = RequestType(alias)

            if alias == RequestType.UNREQUESTED or alias == RequestType.WARN:
                continue
            elif alias == RequestType.REQUESTED and prop in args:
                res[prop] = args[prop]
            elif alias == RequestType.ERROR_IF_PASSED and prop in args:
                raise ValueError(
                    f"{prop} is passed but is not explicitly set as "
                    f"requested or not. In method: {self.name}"
                )
            elif alias in args:
                res[prop] = args[alias]
        return res

    @classmethod
    def deserialize(cls, obj, name):
        """Deserialize an instance from the given object.

        Parameters
        ----------
        obj : dict
            A serialized version of an instance of this object.

        name : str
            The name of the method to which `obj` belongs.

        Returns
        -------
        result : MethodMetdataRequest
            An instance of this class constructed from the input.
        """
        # never change what's passed here.
        requests = deepcopy(obj)
        if requests is None:
            requests = dict()
        result = cls(name=name)
        for prop, alias in requests.items():
            result.add_request(prop=prop, alias=alias)
        return result

    def serialize(self):
        """Serialize the object."""
        result = dict()
        # Then parameters with string aliases
        result.update(
            {
                prop: alias
                for prop, alias in self._requests.items()
                if isinstance(alias, str)
            }
        )
        # And at last the parameters with RequestType routing info
        result.update(
            {
                prop: RequestType(alias).value
                for prop, alias in self._requests.items()
                if not isinstance(alias, str)
            }
        )
        return result

    def __repr__(self):
        return str(self.serialize())

    def __str__(self):
        return str(repr(self))


class MetadataRequest:
    """Contains the metadata request info of a consumer.

    Instances of :class:`MethodMetadataRequest` are used in this class for each
    available method under `metadatarequest.{method}`.

    Consumers-only classes such as simple estimators return a serialized
    version of this class as the output of `get_metadata_routing()`.

    .. versionadded:: 1.1
    """

    def __init__(self):
        for method in METHODS:
            setattr(self, method, MethodMetadataRequest(name=method))

    def _get_param_names(self, method, original_names):
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

        Returns
        -------
        names : set of str
            Returns a set of strings with the names of all parameters.
        """
        return getattr(self, method)._get_param_names(original_names=original_names)

    def _get_params(self, *, method, params):
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
        return getattr(self, method)._get_params(params=params)

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

    def serialize(self):
        """Return dictionary representation of this object."""
        output = {"^type": "request"}
        for method in METHODS:
            output[method] = getattr(self, method).serialize()
        return output

    @classmethod
    def deserialize(cls, obj):
        """Deserialize an instance from the given object.

        Parameters
        ----------
        obj : dict
            A serialized version of an instance of this object.

        Returns
        -------
        result : MetdataRequest
            An instance of this class constructed from the input.
        """
        # never change what's passed here.
        requests = deepcopy(obj)
        result = cls()
        obj_type = requests.pop("^type", None)
        if obj_type != "request":
            raise ValueError(
                "Can only create a router of type 'router', given `_type` is:"
                f" {obj_type}."
            )

        for method, method_requests in requests.items():
            if method not in METHODS:
                raise ValueError(f"{method} is not supported as a method.")
            setattr(
                result,
                method,
                MethodMetadataRequest.deserialize(method_requests, name=method),
            )
        return result

    def __repr__(self):
        return str(self.serialize())

    def __str__(self):
        return str(repr(self))


def metadata_request_factory(obj=None):
    """Get a MetadataRequest instance from the given object.

    This function always returns a copy or an instance constructed from the
    intput, such that changing the output of this function will not change the
    original object.

    .. versionadded:: 1.1

    Parameters
    ----------
    obj : object
        If the object is already a MetadataRequest, return that.
        If the object is an estimator, try to call `get_metadata_request` and get
        an instance from that method.
        If the object is a dict, create a MetadataRequest from that.

    Returns
    -------
    metadata_requests : MetadataRequest
        A ``MetadataRequest`` taken or created from the given object.
    """
    if obj is None:
        return MetadataRequest()

    # doing this instead of a try/except since an AttributeError could be raised
    # for other reasons.
    if hasattr(obj, "get_metadata_routing"):
        return metadata_request_factory(obj.get_metadata_routing())

    if isinstance(obj, (MetadataRouter, MetadataRequest)):
        return deepcopy(obj)

    if isinstance(obj, dict):
        obj_type = obj.get("^type", None)
        if obj_type == "request":
            return MetadataRequest.deserialize(obj)
        elif obj_type == "router":
            return MetadataRouter.deserialize(obj)
        else:
            raise ValueError(f"Cannot understand object type: {obj_type}")

    return MetadataRequest()


class MethodMapping:
    """Stores the mapping from an object's methods to a router's methods.

    This class is primarily used in a ``get_metadata_routing()`` of a router
    object when defining the mapping between a sub-object (a sub-estimator or a
    scorer) to the router's methods. It stores a collection of ``Route``
    namedtuples.

    .. versionadded: 1.1
    """

    def __init__(self):
        self._routes = []

    def __iter__(self):
        for route in self._routes:
            yield (route.method, route.used_in)

    def add(self, *, method, used_in):
        """Add a method mapping.

        Parameters
        ----------
        method : str
            Child object's method name.

        used_in : str
            Parent estimator's method name in which the `"method"` is used.

        Returns
        -------
        self : MethodMapping
            Returns self.
        """
        if method not in METHODS:
            raise ValueError(f"Given method:{method} is not valid.")
        if used_in not in METHODS:
            raise ValueError(f"Given used_in method:{used_in} is not valid.")
        self._routes.append(Route(method=method, used_in=used_in))
        return self

    def serialize(self):
        """Serialize the instance."""
        result = list()
        for route in self._routes:
            result.append({"method": route.method, "used_in": route.used_in})
        return result

    @classmethod
    def deserialize(cls, obj):
        """Deserialize an instance from the given object.

        Parameters
        ----------
        obj : dict
            A serialized version of an instance of this object.

        Returns
        -------
        result : MethodMapping
            An instance of this class constructed from the input.
        """
        # never change what's passed here.
        obj = deepcopy(obj)
        result = cls()
        for route in obj:
            result.add(method=route["method"], used_in=route["used_in"])
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
                routing.add(method=method, used_in=method)
        elif route in METHODS:
            routing.add(method=route, used_in=route)
        else:
            raise ValueError("route should be 'one-to-one' or a single method!")
        return routing

    def __repr__(self):
        return str(self.serialize())


class MetadataRouter:
    """Stores and handles metadata routing for a router object.

    This class is used by router objects to store and handle metadata routing.
    Routing information is stored as a dictionary of the form ``{"object_name":
    RouteMappingPair(method_mapping, routing_info)}``, where ``method_mapping``
    is an instance of :class:`~utils.metadata_requests.MethodMapping` and
    ``routing_info`` is either a
    :class:`~utils.metadata_requests.MetadataRequest` or a
    :class:`~utils.metadata_requests.MetadataRouter` instance.

    .. versionadded: 1.1
    """

    def __init__(self):
        self._route_mappings = dict()

    def _add(self, *, method_mapping, **objs):
        """Add a route mapping w/o name validation.

        Parameters
        ----------
        method_mapping : MethodMapping or str
            The mapping between the child and the parent's methods. If str, the
            output of :func:`~utils.metadata_requests.MethodMapping.from_str`
            is used.

        **objs : dict
            A dictionary of objects from which metadata is extracted by calling
            :func:`~utils.metadata_requests.metadata_request_factory` on them.
        """
        if isinstance(method_mapping, str):
            method_mapping = MethodMapping.from_str(method_mapping)
        for name, obj in objs.items():
            self._route_mappings[name] = RouteMappingPair(
                mapping=method_mapping, routing=metadata_request_factory(obj)
            )

    def add_self(self, obj):
        """Add `self` to the routing map.

        A "one-to-one" method map is used for method mapping.

        Parameters
        ----------
        obj : object
            This is typically `self` in a `get_metadata_routing()`.
        """
        self._add(method_mapping="one-to-one", me=obj._get_metadata_request())
        return self

    def add(self, *, method_mapping, **objs):
        """Add named objects with their corresponding method mapping.

        The objects cannot be called `"me"` since that's reserved for `self`.

        Parameters
        ----------
        method_mapping : MethodMapping or str
            The mapping between the child and the parent's methods. If str, the
            output of :func:`~utils.metadata_requests.MethodMapping.from_str`
            is used.

        **objs : dict
            A dictionary of objects from which metadata is extracted by calling
            :func:`~utils.metadata_requests.metadata_request_factory` on them.
        """
        if "me" in objs.keys():
            raise ValueError("'me' is reserved for `self`! Use a different keyword")
        self._add(method_mapping=method_mapping, **objs)
        return self

    def _get_self(self):
        """Return the self metadata request if available.

        If a router is also a consumer, this information is stored under the
        reserved key `"me"`.
        """
        if "me" in self._route_mappings:
            return self._route_mappings["me"].routing
        else:
            return MetadataRequest()

    def _get_param_names(self, *, method, original_names):
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

        Returns
        -------
        names : set of str
            Returns a set of strings with the names of all parameters.
        """
        res = set()
        for name, route_mapping in self._route_mappings.items():
            # if original names are required, that only applies to "me"
            orig_names = original_names and name == "me"
            for orig_method, used_in in route_mapping.mapping:
                if used_in == method:
                    res = res.union(
                        route_mapping.routing._get_param_names(
                            method=orig_method, original_names=orig_names
                        )
                    )
        return set(sorted(res))

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
        param_names = self._get_param_names(method=method, original_names=False)
        return Bunch(
            **{key: value for key, value in params.items() if key in param_names}
        )

    def get_params(self, *, method, params):
        """Return the input parameters requested by a child objects.

        The output of this method is a bunch, which includes the inputs for all
        methods of each child object that are used in the router's `"method"`.

        If the router is also a consumer, it also checks for warnings of
        `self`'s requested metadata.

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
            A :class:`~utils.Bunch` of the form
            ``{"object_name": {"method_name": {prop: value}}}`` which can be
            used to pass the required metadata to corresponding methods or
            corresponding child objects.
        """
        res = Bunch()
        for name, route_mapping in self._route_mappings.items():
            router, mapping = route_mapping.routing, route_mapping.mapping
            if name == "me":
                # "me" is reserved for `self` and routing is not handled by
                # `self`, it's handled by meta-estimators. Therefore we only
                # check for warnings here.
                router._check_warnings(params=params, method=method)
                continue

            res[name] = Bunch()
            for orig_method, used_in in mapping:
                if used_in == method:
                    if isinstance(router, MetadataRequest):
                        res[name][orig_method] = router._get_params(
                            params=params, method=orig_method
                        )
                    else:
                        res[name][orig_method] = router._get_squashed_params(
                            params=params, method=orig_method
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
            routed.

        params : dict
            A dictionary of provided metadata.
        """
        param_names = self._get_param_names(method=method, original_names=True)
        self_params = self._get_self()._get_param_names(
            method=method, original_names=True
        )
        extra_keys = set(sorted(set(params.keys()) - param_names - self_params))
        if extra_keys:
            raise ValueError(
                "These passed parameters are not understood or requested by any object:"
                f" {extra_keys}"
            )

    def serialize(self):
        """Serialize the instance."""
        res = {"^type": "router"}
        for name, route_mapping in self._route_mappings.items():
            res[name] = dict()
            res[name]["mapping"] = route_mapping.mapping.serialize()
            res[name]["routing"] = route_mapping.routing.serialize()

        return res

    @classmethod
    def deserialize(cls, obj):
        """Deserialize an instance from the given object.

        Parameters
        ----------
        obj : dict
            A serialized version of an instance of this object.

        Returns
        -------
        result : MetadataRouter
            An instance of this class constructed from the input.
        """
        # never change what's passed here.
        obj = deepcopy(obj)
        obj_type = obj.pop("^type", None)
        if obj_type != "router":
            raise ValueError(
                "Can only create a router of type 'router', given `_type` is:"
                f" {obj_type}."
            )

        res = cls()
        for name, route_mapping in obj.items():
            res._add(
                method_mapping=MethodMapping.deserialize(route_mapping["mapping"]),
                **{name: metadata_request_factory(route_mapping["routing"])},
            )
        return res

    def __iter__(self):
        for name, route_mapping in self._route_mappings.items():
            yield (name, route_mapping)


class RequestMethod:
    """
    A descriptor for request methods.

    .. versionadded:: 1.1

    Parameters
    ----------
    name : str
        The name of the method for which the request function should be
        created, e.g. ``"fit"`` would create a ``fit_requests`` function.

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
                raise TypeError(f"Unexpected args: {set(kw) - set(self.keys)}")

            requests = metadata_request_factory(instance._get_metadata_request())

            try:
                method_metadata_request = getattr(requests, self.name)
            except AttributeError:
                raise ValueError(f"{self.name} is not a supported method.")

            for prop, alias in kw.items():
                if alias is not UNCHANGED:
                    method_metadata_request.add_request(prop=prop, alias=alias)
            instance._metadata_request = requests.serialize()

            return instance

        # Now we set the relevant attributes of the function so that it seems
        # like a normal method to the end user, with known expected arguments.
        func.__name__ = f"{self.name}_requests"
        params = [
            inspect.Parameter(
                name="self",
                kind=inspect.Parameter.POSITIONAL_OR_KEYWORD,
                annotation=type(instance),
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
            return_annotation=type(instance),
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
        """Set the ``{method}_requests`` methods.

        This uses PEP-487 [1]_ to set the ``{method}_requests`` methods. It
        looks for the information available in the set default values which are
        set using ``__metadata_request__*`` class attributes.

        References
        ----------
        .. [1] https://www.python.org/dev/peps/pep-0487
        """
        try:
            requests = cls._get_default_requests().serialize()
        except Exception:
            # if there are any issues in the default values, it will be raised
            # when ``get_metadata_routing`` is called. Here we are going to
            # ignore all the issues such as bad defaults etc.
            super().__init_subclass__(**kwargs)
            return

        for request_method, request_keys in requests.items():
            # ignore everything which is not a valid method, including
            # serialized metadata such as ^type
            if request_method not in METHODS:
                continue
            # set ``{method}_requests``` methods
            if not len(request_keys):
                continue
            setattr(
                cls,
                f"{request_method}_requests",
                RequestMethod(request_method, sorted(request_keys)),
            )
        super().__init_subclass__(**kwargs)

    @classmethod
    def _get_default_requests(cls):
        """Collect default request values.

        This method combines the information present in ``metadata_request__*``
        class attributes.
        """

        requests = MetadataRequest()

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
        request : dict
            A dict of dict of str->value. The key to the first dict is the name
            of the method, and the key to the second dict is the name of the
            argument requested by the method.
        """
        if hasattr(self, "_metadata_request"):
            requests = metadata_request_factory(self._metadata_request)
        else:
            requests = self._get_default_requests()

        return requests.serialize()

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


def process_routing(func):
    """Decorator to process input params and handle routing.

    This decorator wraps around a router estimator's method, such as ``fit``,
    and handles routing and and validation of the metadata passed to it.

    .. versionadded: 1.1

    Parameters
    ----------
    func : method
        A router estimator's method such as ``fit``.

    Returns
    -------
    decorator : method
        The estimator's method, after validation and routing of the metadata.
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        sig = inspect.signature(func)
        # Create the argument binding so we can determine what
        # parameters are given what values
        argument_binding = sig.bind(*args, **kwargs)
        argument_map = argument_binding.arguments

        # find parameter names which should be passed explicitly, and the kwarg param
        # which can accept routed parameters.
        func_params = list(sig.parameters.items())
        params_to_pass = []
        kwarg_param = None
        for pname, param in func_params:
            if param.kind == param.VAR_KEYWORD:
                kwarg_param = pname
                continue
            params_to_pass.append(pname)

        if kwarg_param is None:
            raise TypeError("Method must accept kwargs for metadata routing to work.")

        obj = argument_map["self"]
        if not hasattr(obj, "get_metadata_routing"):
            raise AttributeError(
                f"This {repr(obj.__class__.__name__)} has not implemented the routing"
                " method `get_metadata_routing`."
            )
        method_name = func.__name__
        if method_name not in METHODS:
            raise TypeError(
                f"Can only route and process input on these methods: {METHODS}, "
                f"while the current method is: {method_name}."
            )

        considered_params = {
            key: value
            for key, value in argument_map.items()
            if key not in ("X", "y", "Y", "self")
        }
        # WWe pop the kwarg parameter from the params above, and then update it
        # with the remaining values. This resembles the code snippets such as:
        # if sample_weight is not None:
        #     fit_params["sample_weight"] = sample_weight
        if kwarg_param in considered_params:
            routed_params = considered_params.pop(kwarg_param)
            routed_params.update(considered_params)
        else:
            routed_params = considered_params

        request_routing = metadata_request_factory(obj)
        request_routing.validate_metadata(params=routed_params, method=method_name)
        params = request_routing.get_params(params=routed_params, method=method_name)

        # create a dict of other parameters and then add 'params' to them
        explicit_params = {
            key: value for key, value in argument_map.items() if key in params_to_pass
        }
        explicit_params["params"] = params
        return func(**explicit_params)

    return wrapper
