# from copy import deepcopy
from enum import Enum
from ..externals._sentinels import sentinel  # type: ignore # mypy error!!!


class RequestType(Enum):
    UNREQUESTED = False
    REQUESTED = True
    ERROR_IF_PASSED = None


UNCHANGED = sentinel("UNCHANGED")

METHODS = [
    "fit",
    "partial_fit",
    "predict",
    "score",
    "split",
    "transform",
    "inverse_transform",
]


class MethodMetadataRequest:
    """Contains the metadata request info for a single method."""

    def __init__(self):
        self.requests = dict()

    def add_request(
        self,
        *,
        prop,
        alias,
        allow_aliasing=True,
        overwrite=False,
        expected_metadata=None,
    ):
        """Add request info for a prop.

        Parameters
        ----------
        prop : str
            The property for which a request is set.

        alias : str or {True, False, None}
            The alias which is routed to `prop`
            True: requested
            False: not requested
            None: error if passed

        allow_aliasing : bool, default=True
            If False, alias should be the same as prop if it's a string.

        overwrite : bool or str, default=False
            - True: ``alias`` replaces the existing routing.
            - False: a ``ValueError`` is raised if the given value conflicts
              with an existing one.
            - "on-default": only overwrite the alias is it is
              RequestType.ERROR_IF_PASSED

        expected_metadata : str, default=None
            If provided, all props should be the same as this value. It used to
            handle default values.
        """
        if overwrite not in {True, False, "on-default"}:
            raise ValueError(
                "overwrite can only be one of {True, False, 'on-default'}."
            )
        if expected_metadata is not None and expected_metadata != prop:
            raise ValueError(
                f"Expected all metadata to be called {expected_metadata} but "
                f"{prop} was passed."
            )
        if not allow_aliasing and isinstance(alias, str) and prop != alias:
            raise ValueError(
                "Aliasing is not allowed, prop and alias should "
                "be the same strings if alias is a string."
            )

        if not isinstance(alias, str):
            try:
                RequestType(alias)
            except ValueError:
                raise ValueError(
                    "alias should be either a string or one of "
                    "{None, True, False}, or a RequestType."
                )

        if alias == prop:
            alias = RequestType.REQUESTED
        if not isinstance(alias, str):
            alias = RequestType(alias)

        if prop not in self.requests or overwrite:
            self.requests[prop] = alias
        elif (
            overwrite == "on-default"
            and not isinstance(self.requests[prop], str)
            and RequestType(self.requests[prop]) == RequestType.ERROR_IF_PASSED
        ):
            self.requests[prop] = alias
        elif self.requests[prop] != alias:
            raise ValueError(
                f"{prop} is already requested as {self.requests[prop]}, "
                f"which is not the same as the one given: {alias}."
            )

    def merge_method_request(self, other, overwrite=False, expected_metadata=None):
        """Merge the metadata request info of two methods.

        The methods can be the same, or different. For example, merging
        fit and score info of the same object, or merging fit request info
        from two different sub estimators.

        Parameters
        ----------
        other : MethodMetadataRequest
            The other object to be merged with this instance.

        overwrite : bool or str, default=False
            - True: ``alias`` replaces the existing routing.
            - False: a ``ValueError`` is raised if the given value conflicts
              with an existing one.
            - "on-default": only overwrite the alias is it is
              RequestType.ERROR_IF_PASSED

        expected_metadata : str, default=None
            If provided, all props should be the same as this value. It used to
            handle default values.
        """
        if not isinstance(other, MethodMetadataRequest):
            raise ValueError("Can only add another MethodMetadataRequest.")
        for prop, alias in other.requests.items():
            self.add_request(
                prop=prop,
                alias=alias,
                overwrite=overwrite,
                expected_metadata=expected_metadata,
            )

    def validate_metadata(self, ignore_extras=False, **kwargs):
        """Validate the given arguments against the requested ones.

        Parameters
        ----------
        ignore_extras : bool, default=False
            If ``True``, no error is raised if extra unknown args are passed.

        kwargs : dict
            Provided metadata.

        Returns
        -------
        None
        """
        args = {arg for arg, value in kwargs.items() if value is not None}
        if not ignore_extras and args - set(self.requests.keys()):
            raise ValueError(
                "Metadata passed which is not understood: "
                f"{args - set(self.requests.keys())}"
            )

        for prop, alias in self.requests.items():
            if not isinstance(alias, str):
                alias = RequestType(alias)
            if alias == RequestType.UNREQUESTED:
                continue
            elif alias == RequestType.REQUESTED or isinstance(alias, str):
                # we ignore what the given alias here is, since aliases are
                # checked at the parent meta-estimator level, and the child
                # still expects the original names for the metadata.
                # If a metadata is requested but not passed, no error is raised
                continue
            elif alias == RequestType.ERROR_IF_PASSED:
                if prop in args:
                    raise ValueError(
                        f"{prop} is passed but is not explicitly set as "
                        "requested or not."
                    )

    def get_method_input(self, ignore_extras=False, **kwargs):
        """Return the input parameters requested by the method.

        The output of this method can be used directly as the input to the
        corresponding method as extra props.

        Parameters
        ----------
        ignore_extras : bool, default=False
            If ``True``, no error is raised if extra unknown args are passed.

        kwargs : dict
            A dictionary of provided metadata.

        Returns
        -------
        kwargs : dict
            A dictionary of {prop: value} which can be given to the
            corresponding method.
        """
        args = {arg: value for arg, value in kwargs.items() if value is not None}
        res = dict()
        for prop, alias in self.requests.items():
            if not isinstance(alias, str):
                alias = RequestType(alias)

            if alias == RequestType.UNREQUESTED:
                continue
            elif alias in {RequestType.REQUESTED, RequestType.ERROR_IF_PASSED}:
                # if alias is ERROR_IS_PASSED, putting it in res makes sure the
                # call to validate will raise the appropriate error.
                # Also, if prop is not present in args, putting None raises the
                # right error in validate_metadata.
                res[prop] = args.get(prop, None)
            else:
                res[prop] = args.get(alias, None)
        self.validate_metadata(ignore_extras=ignore_extras, **res)
        return res

    def masked(self):
        """Return a masked version of the requests.

        A masked version is one which converts a ``{'prop': 'alias'}`` to
        ``{'alias': True}``. This is desired in meta-estimators passing
        requests to their parent estimators.
        """
        res = MethodMetadataRequest()
        for prop, alias in self.requests.items():
            if isinstance(alias, str):
                res.add_request(
                    prop=alias,
                    alias=alias,
                    allow_aliasing=False,
                    overwrite=False,
                )
            else:
                res.add_request(
                    prop=prop,
                    alias=alias,
                    allow_aliasing=False,
                    overwrite=False,
                )
        return res

    def to_dict(self):
        """Dictionary representation of this object."""
        return self.requests

    @classmethod
    def from_dict(cls, requests, allow_aliasing=True):
        """Construct a MethodMetadataRequest from a given dictionary.

        Parameters
        ----------
        requests : dict
            A dictionary representing the requests.

        allow_aliasing : bool, default=True
            If false, only aliases with the same name as the parameter are
            allowed. This is useful when handling the default values.

        Returns
        -------
        requests: MethodMetadataRequest
            A :class:`MethodMetadataRequest` object.
        """
        if requests is None:
            requests = dict()
        elif isinstance(requests, str):
            requests = {requests: requests}
        elif isinstance(requests, list):
            requests = {r: r for r in requests}
        result = cls()
        for prop, alias in requests.items():
            result.add_request(prop=prop, alias=alias, allow_aliasing=allow_aliasing)
        return result

    def __repr__(self):
        return self.to_dict()

    def __str__(self):
        return str(self.to_dict())


class MetadataRequest:
    """Contains the metadata request info of an object."""

    def __init__(self, obj=None):
        for method in METHODS:
            setattr(self, method, MethodMetadataRequest())

        if obj is None:
            return
        elif not isinstance(obj, dict):
            raise ValueError(
                "Can only construct an instance from a dict. Please call "
                "metadata_request_factory for other types of input."
            )

        for method, requests in obj.items():
            try:
                mmr = getattr(self, method)
            except AttributeError:
                raise ValueError(f"{method} is not supported as a method.")
            if isinstance(requests, str):
                requests = {requests: None}
            for prop, alias in requests.items():
                mmr.add_request(prop=prop, alias=alias)

    def add_requests(
        self,
        obj,
        mapping="one-to-one",
        overwrite=False,
        expected_metadata=None,
    ):
        """Add request info from the given object with the desired mapping.

        Parameters
        ----------
        obj : object
            An object from which a MetadataRequest can be constructed.

        mapping : dict or str, default="one-to-one"
            The mapping between the ``obj``'s methods and this object's
            methods. If ``"one-to-one"`` all methods' requests from ``obj`` are
            merged into this instance's methods.

        overwrite : bool or str, default=False
            - True: ``alias`` replaces the existing routing.
            - False: a ``ValueError`` is raised if the given value conflicts
              with an existing one.
            - "on-default": only overwrite the alias is it is
              RequestType.ERROR_IF_PASSED

        expected_metadata : str, default=None
            If provided, all props should be the same as this value. It used to
            handle default values.
        """
        if not isinstance(mapping, dict) and mapping != "one-to-one":
            raise ValueError(
                "mapping can only be a dict or the literal 'one-to-one'. "
                f"Given value: {mapping}"
            )
        if mapping == "one-to-one":
            mapping = {method: method for method in METHODS}
        other = metadata_request_factory(obj)
        for destination, source in mapping.items():
            my_method = getattr(self, destination)
            other_method = getattr(other, source)
            my_method.merge_method_request(
                other_method,
                overwrite=overwrite,
                expected_metadata=expected_metadata,
            )

    def masked(self):
        """Return a masked version of the requests.

        A masked version is one which converts a ``{'prop': 'alias'}`` to
        ``{'alias': True}``. This is desired in meta-estimators passing
        requests to their parent estimators.
        """
        res = MetadataRequest()
        for method in METHODS:
            setattr(res, method, getattr(self, method).masked())
        return res

    def to_dict(self):
        """Dictionary representation of this object."""
        output = dict()
        for method in METHODS:
            output[method] = getattr(self, method).to_dict()
        return output

    def __repr__(self):
        return self.to_dict()

    def __str__(self):
        return str(self.to_dict())


def metadata_request_factory(obj=None):
    """Get a MetadataRequest instance from the given object.

    If the object is already a MetadataRequest, return that.
    If the object is an estimator, try to call `get_metadata_request` and get
    an instance from that method.
    If the object is a dict, create a MetadataRequest from that.
    """
    if obj is None:
        return MetadataRequest()

    if isinstance(obj, MetadataRequest):
        return obj

    if isinstance(obj, dict):
        return MetadataRequest(obj)

    try:
        return obj.get_metadata_request(output="MetadataRequest")
    except TypeError:
        # The method exists, but doesn't accept `output`
        return MetadataRequest(obj.get_metadata_request())
    except AttributeError:
        # The object doesn't have a `get_metadata_request` method.
        return MetadataRequest()


class MetadataRouter:
    def __init__(self):
        self.requests = MetadataRequest()

    def add(self, *obj, mapping="one-to-one", overwrite=False, mask=False):
        """Add a set of requests to the existing ones.

        Parameters
        ----------
        *obj : objects
            A set of objects from which the requests are extracted. Passed as
            arguments to this method.

        mapping : dict or str, default="one-to-one"
            The mapping between the ``obj``'s methods and this object's
            methods. If ``"one-to-one"`` all methods' requests from ``obj`` are
            merged into this instance's methods.

        overwrite : bool or str, default=False
            - True: ``alias`` replaces the existing routing.
            - False: a ``ValueError`` is raised if the given value conflicts
              with an existing one.
            - "on-default": only overwrite the alias is it is
              RequestType.ERROR_IF_PASSED

        mask : bool, default=False
            If the requested metadata should be masked by the alias. If
            ``True``, then a request of the form
            ``{'sample_weight' : 'my_weight'}`` is converted to
            ``{'my_weight': my_weight'}``. This is required for meta-estimators
            which should expose the requested parameters and not the ones
            expected by the objects' methods.
        """
        for x in obj:
            if mask:
                x = metadata_request_factory(x).masked()
            self.requests.add_requests(x, mapping=mapping, overwrite=overwrite)
        return self

    def get_metadata_request(self, output="dict"):
        """Get requested data properties.

        Parameters
        ----------
        output : {"dict", "MetadataRequest}
            Whether the output should be a MetadataRequest instance, or a dict
            representing that instance.

        Returns
        -------
        request : MetadataRequest, or dict
            If dict, it will be a deserialized version of the underlying
            MetadataRequest object: dict of dict of str->value. The key to the
            first dict is the name of the method, and the key to the second
            dict is the name of the argument requested by the method.
        """
        if output == "dict":
            return self.requests.to_dict()
        elif output == "MetadataRequest":
            return self.requests
        else:
            raise ValueError("output can be one of {'dict', 'MetadataRequest'}")
