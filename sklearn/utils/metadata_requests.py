# from copy import deepcopy


class Target:
    def __init__(self, *, owner, name, is_default):
        self.owner = owner
        self.name = name
        self.is_default = is_default


class MethodMetadataRequest:
    def __init__(self):
        self.requests = dict()

    def add_request(self, alias, targets):
        if alias not in self.requests:
            self.requests[alias] = set()
        self.requests[alias] = self.requests[alias].union(targets)

    def add_method_request(self, other):
        if not isinstance(other, MethodMetadataRequest):
            raise ValueError("Can only add another MethodMetadataRequest.")
        for alias, targets in other.requests.items():
            self.add_request(alias, targets)


class MetadataRequest:
    def __init__(self, obj=None):
        metadata_request = obj.get_metadata_request()
        if isinstance(metadata_request, MetadataRequest):
            return metadata_request

        self.fit = MethodMetadataRequest()
        self.predict = MethodMetadataRequest()
        self.score = MethodMetadataRequest()
        self.split = MethodMetadataRequest()
        self.transform = MethodMetadataRequest()

        if obj is None:
            return self

        for method, requests in obj.get_metadata_request():
            mmr = getattr(self, method)
            for alias, target in requests.items():
                target = Target(
                    owner=obj,
                    name=target if target is not None else alias,
                    is_default=target is None,
                )
                mmr.add_request(alias, {target})
        return self

    def add_requests(obj, mapping=None):
        if mapping is None:
            mapping = {
                "fit": "fit",
                "predict": "predict",
                "score": "score",
                "split": "split",
                "transform": "transform",
            }
