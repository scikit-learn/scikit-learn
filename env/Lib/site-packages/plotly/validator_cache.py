import importlib
from _plotly_utils.basevalidators import LiteralValidator


class ValidatorCache(object):
    _cache = {}

    @staticmethod
    def get_validator(parent_path, prop_name):

        key = (parent_path, prop_name)
        if key not in ValidatorCache._cache:

            if "." not in parent_path and prop_name == "type":
                # Special case for .type property of traces
                validator = LiteralValidator("type", parent_path, parent_path)
            else:
                lookup_name = None
                if parent_path == "layout":
                    from .graph_objects import Layout

                    match = Layout._subplotid_prop_re.match(prop_name)
                    if match:
                        lookup_name = match.group(1)

                lookup_name = lookup_name or prop_name
                class_name = lookup_name.title() + "Validator"
                validator = getattr(
                    importlib.import_module("plotly.validators." + parent_path),
                    class_name,
                )(plotly_name=prop_name)
            ValidatorCache._cache[key] = validator

        return ValidatorCache._cache[key]
