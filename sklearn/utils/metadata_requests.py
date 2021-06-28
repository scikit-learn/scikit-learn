# from copy import deepcopy
import inspect
from enum import Enum
from collections import defaultdict
from typing import Union, Optional
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
    """Contains the metadata request info for a single method.

    .. versionadded:: 1.1

    Parameters
    ----------
    name : str
        The name of the method to which these requests belong.
    """

    def __init__(self, name):
        self.requests = dict()
        self.name = name

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

    def validate_metadata(self, ignore_extras=False, kwargs=None):
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
        kwargs = {} if kwargs is None else kwargs
        args = {arg for arg, value in kwargs.items() if value is not None}
        if not ignore_extras and args - set(self.requests.keys()):
            raise ValueError(
                "Metadata passed which is not understood: "
                f"{args - set(self.requests.keys())}. In method: {self.name}"
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
                        f"requested or not. In method: {self.name}"
                    )

    def get_method_input(self, ignore_extras=False, kwargs=None):
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
        kwargs = {} if kwargs is None else kwargs
        args = {arg: value for arg, value in kwargs.items() if value is not None}
        res = dict()
        for prop, alias in self.requests.items():
            if not isinstance(alias, str):
                alias = RequestType(alias)

            if alias == RequestType.UNREQUESTED:
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
        self.validate_metadata(ignore_extras=ignore_extras, kwargs=res)
        return res

    def masked(self):
        """Return a masked version of the requests.

        Returns
        -------
        masked : MethodMetadataRequest
            A masked version is one which converts a ``{'prop': 'alias'}`` to
            ``{'alias': True}``. This is desired in meta-estimators passing
            requests to their parent estimators.
        """
        res = MethodMetadataRequest(name=self.name)
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

    @classmethod
    def from_dict(cls, requests, name, allow_aliasing=True):
        """Construct a MethodMetadataRequest from a given dictionary.

        Parameters
        ----------
        requests : dict
            A dictionary representing the requests.

        name : str
            The name of the method to which these requests belong.

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
            requests = {requests: RequestType.ERROR_IF_PASSED}
        elif isinstance(requests, list):
            requests = {r: RequestType.ERROR_IF_PASSED for r in requests}
        result = cls(name=name)
        for prop, alias in requests.items():
            result.add_request(prop=prop, alias=alias, allow_aliasing=allow_aliasing)
        return result

    def __repr__(self):
        return self.requests

    def __str__(self):
        return str(self.requests)


class MetadataRequest:
    """Contains the metadata request info of an object.

    .. versionadded:: 1.1

    Parameters
    ----------
    requests : dict of dict of {str: str}, default=None
        A dictionary where the keys are the names of the methods, and the values are
        a dictionary of the form ``{"required_metadata": "provided_metadata"}``.
        ``"provided_metadata"`` can also be a ``RequestType`` or {True, False, None}.

    default : RequestType, True, False, None, or str, \
        default=RequestType.ERROR_IF_PASSED
        The default value to be used if parameters are provided as a string instead of
        the usual second layer dict.
    """

    def __init__(self, requests=None, default=RequestType.ERROR_IF_PASSED):
        for method in METHODS:
            setattr(self, method, MethodMetadataRequest(name=method))

        if requests is None:
            return
        elif not isinstance(requests, dict):
            raise ValueError(
                "Can only construct an instance from a dict. Please call "
                "metadata_request_factory for other types of input."
            )

        for method, method_requests in requests.items():
            method_requests = {} if method_requests is None else method_requests
            try:
                mmr = getattr(self, method)
            except AttributeError:
                raise ValueError(f"{method} is not supported as a method.")
            if isinstance(method_requests, str):
                method_requests = {method_requests: default}
            elif isinstance(method_requests, (list, set)):
                method_requests = {m: default for m in method_requests}
            for prop, alias in method_requests.items():
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
            merged into this instance's methods. If a dict, the mapping is of
            the form ``{"destination_method": "source_method"}``.

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
        """Return dictionary representation of this object."""
        output = dict()
        for method in METHODS:
            output[method] = getattr(self, method).requests
        return output

    def __repr__(self):
        return self.to_dict()

    def __str__(self):
        return str(self.to_dict())


def metadata_request_factory(obj=None):
    """Get a MetadataRequest instance from the given object.

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

    if isinstance(obj, MetadataRequest):
        return obj

    if isinstance(obj, dict):
        return MetadataRequest(obj)

    try:
        return MetadataRequest(obj.get_metadata_request())
    except AttributeError:
        # The object doesn't have a `get_metadata_request` method.
        return MetadataRequest()


class MetadataRouter:
    """Route the metadata to child objects.

    .. versionadded:: 1.1
    """

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
            The mapping between the ``obj``'s methods and this routing object's
            methods. If ``"one-to-one"`` all methods' requests from ``obj`` are
            merged into this instance's methods. If a dict, the mapping is of
            the form ``{"destination_method": "source_method"}``.

        overwrite : bool or str, default=False
            - True: ``alias`` replaces the existing routing.
            - False: a ``ValueError`` is raised if the given value conflicts
              with an existing one.
            - "on-default": only overwrite the alias if it is
              RequestType.ERROR_IF_PASSED

        mask : bool, default=False
            If the requested metadata should be masked by the alias. If
            ``True``, then a request of the form
            ``{'sample_weight' : 'my_weight'}`` is converted to
            ``{'my_weight': 'my_weight'}``. This is required for meta-estimators
            which should expose the requested parameters and not the ones
            expected by the objects' methods.
        """
        for x in obj:
            if mask:
                x = metadata_request_factory(x).masked()
            self.requests.add_requests(x, mapping=mapping, overwrite=overwrite)
        return self

    def get_metadata_request(self):
        """Get requested data properties.

        Returns
        -------
        request : dict
            A dict of dict of str->value. The key to the first dict is the name
            of the method, and the key to the second dict is the name of the
            argument requested by the method.
        """
        return self.requests.to_dict()


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

            requests = metadata_request_factory(instance)

            try:
                method_metadata_request = getattr(requests, self.name)
            except AttributeError:
                raise ValueError(f"{self.name} is not a supported method.")

            for prop, alias in kw.items():
                if alias is not UNCHANGED:
                    method_metadata_request.add_request(
                        prop=prop, alias=alias, allow_aliasing=True, overwrite=True
                    )
            instance._metadata_request = requests.to_dict()

            return instance

        # Now we set the relevant attributes of the function so that it seems
        # like a normal method to the end user, with known expected arguments.
        func.__name__ = f"{self.name}_requests"
        func.__signature__ = inspect.Signature(
            [
                inspect.Parameter(
                    k,
                    inspect.Parameter.KEYWORD_ONLY,
                    default=UNCHANGED,
                    annotation=Optional[Union[RequestType, str]],
                )
                for k in self.keys
            ],
            return_annotation=type(instance),
        )
        doc = """Request metadata passed to the ``{method}`` method.

            Parameters
            ----------
            """.format(
            method=self.name
        )
        for metadata in self.keys:
            doc += """{metadata} : RequestType, str, True, False, or None, \
                    default=UNCHANGED
                Whether {metadata} should be passed to {method} by meta-estimators or
                not, and if yes, should it have an alias.

                - True or RequestType.REQUESTED: {metadata} is requested, and passed to
                {method} if provided.
                - False or RequestType.UNREQUESTED: {metadata} is not requested and the
                meta-estimator will not pass it to {method}.
                - None or RequestType.ERROR_IF_PASSED: {metadata} is not requested, and
                the meta-estimator will raise an error if the user provides {metadata}
                - str: {metadata} should be passed to the meta-estimator with this given
                alias instead of the original name.

                """.format(
                metadata=metadata, method=self.name
            )
        doc += """
            Returns
            -------
            self : object
                Returns the object itself.
            """
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
        set using ``_metadata_request__*`` class attributes.

        References
        ----------
        .. [1] https://www.python.org/dev/peps/pep-0487
        """
        try:
            requests = cls._get_default_requests().to_dict()
        except Exception:
            # if there are any issues in the default values, it will be raised
            # when ``get_metadata_request`` is called. Here we are going to
            # ignore all the issues such as bad defaults etc.`
            super().__init_subclass__(**kwargs)
            return

        for request_method, request_keys in requests.items():
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
        defaults = defaultdict()
        for klass in reversed(inspect.getmro(cls)):
            klass_defaults = {
                attr: value
                for attr, value in vars(klass).items()
                if attr.startswith("_metadata_request__")
            }
            defaults.update(klass_defaults)
        defaults = dict(sorted(defaults.items()))
        for attr, value in defaults.items():
            requests.add_requests(
                value, expected_metadata="__".join(attr.split("__")[1:])
            )
        return requests

    def get_metadata_request(self):
        """Get requested data properties.

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

        return requests.to_dict()


class SampleWeightConsumer:
    """Mixin class to add ``sample_weight`` request to ``fit`` and ``score``.

    .. versionadded:: 1.1
    """

    _metadata_request__sample_weight = {
        "fit": "sample_weight",
        "score": "sample_weight",
    }
