# from copy import deepcopy
from enum import Enum
from ..externals._sentinels import sentinel # type: ignore # mypy error!!!


class RequestType(Enum):
    UNREQUESTED = False
    REQUESTED = True
    ERROR_IF_PASSED = None


UNCHANGED = sentinel("UNCHANGED")


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
        overwrite=True,
        expected_metadata=None,
    ):
        """Add request info for a prop.

        Parameters
        ----------
        prop: str
            The property for which a request is set.

        alias: str or {True, False, None}
            The alias which is routed to `prop`

        allow_aliasing: bool, default=True
            If False, alias should be the same as prop if it's a string.

        overwrite: bool, default=True
            True: ``alias`` replaces the existing routing.
            False: a ``ValueError`` is raised if the given value conflicts with
            an existing one.

        expected_metadata: str, default=None
            If provided, all props should be the same as this value. It used to
            handle default values.
        """
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

        if alias not in {None, True, False} and not isinstance(alias, str):
            raise ValueError(
                "alias should be either a string or one of "
                "None, True, or False."
            )

        if alias is True:
            alias = prop

        if prop not in self.requests or overwrite:
            self.requests[prop] = alias
        elif self.requests[prop] != alias:
            raise ValueError(
                f"{prop} is already assigned "
                f"an alias as f{self.requests[prop]}, which is "
                f"not the same as the one given: {alias}."
            )

    def merge_method_request(
        self, other, overwrite=False, expected_metadata=None
    ):
        """Merge the metadata request info of two methods.

        The methods can be the same, or different. For example, merging
        fit and score info of the same object, or merging fit request info
        from two different sub estimators.

        Parameters
        ----------
        other: MethodMetadataRequest
            The other object to be merged with this instance.

        overwrite: bool, default=False
            True: ``alias`` replaces the existing routing.
            False: a ``ValueError`` is raised if the given value conflicts with
            an existing one.

        expected_metadata: str, default=None
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

    def to_dict(self):
        """Dictionary representation of this object."""
        return self.requests

    @classmethod
    def from_dict(cls, requests, allow_aliasing=True):
        """Construct a MethodMetadataRequest from a given dictionary.

        Parameters
        ----------
        requests: dict
            A dictionary representing the requests.
        allow_aliasing: bool, default=True
            If false, only aliases with the same name as the parameter are
            allowed. This is useful when handling the default values.

        Returns
        -------
        requests: MethodMetadataRequest
            A :class:`MethodMetadataRequest` object.
        """
        result = cls()
        for prop, alias in requests:
            result.add_request(
                prop=prop, alias=alias, allow_aliasing=allow_aliasing
            )

    def __repr__(self):
        return self.to_dict()

    def __str__(self):
        return str(self.to_dict())


class MetadataRequest:
    """Contains the metadata request info of an object."""

    def __init__(self, obj=None):
        self.fit = MethodMetadataRequest()
        self.partial_fit = MethodMetadataRequest()
        self.predict = MethodMetadataRequest()
        self.score = MethodMetadataRequest()
        self.split = MethodMetadataRequest()
        self.transform = MethodMetadataRequest()
        self.inverse_transform = MethodMetadataRequest()

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
        self, obj, mapping=None, overwrite=False, expected_metadata=None
    ):
        """Add request info from the given object with the desired mapping.

        Parameters
        ----------
        obj: object
            An object from which a MetadataRequest can be constructed.

        mapping: dict, default=None
            The mapping between the ``obj``'s methods and this object's
            methods. By default all methods' requests from ``obj`` are merged
            into this instance's methods.

        overwrite: bool, default=False
            True: ``alias`` replaces the existing routing.
            False: a ``ValueError`` is raised if the given value conflicts with
            an existing one.

        expected_metadata: str, default=None
            If provided, all props should be the same as this value. It used to
            handle default values.
        """
        if mapping is None:
            mapping = {
                "fit": "fit",
                "partial_fit": "partial_fit",
                "predict": "predict",
                "score": "score",
                "split": "split",
                "transform": "transform",
                "inverse_transform": "inverse_transform",
            }
        other = MetadataRequest(obj)
        for source, destination in mapping.items():
            my_method = getattr(self, destination)
            other_method = getattr(other, source)
            my_method.merge_method_request(
                other_method,
                overwrite=overwrite,
                expected_metadata=expected_metadata,
            )

    def to_dict(self):
        """Dictionary representation of this object."""
        output = dict()
        output["fit"] = self.fit.to_dict()
        output["partial_fit"] = self.partial_fit.to_dict()
        output["predict"] = self.predict.to_dict()
        output["score"] = self.score.to_dict()
        output["split"] = self.split.to_dict()
        output["transform"] = self.transform.to_dict()
        output["inverse_transform"] = self.inverse_transform.to_dict()
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
