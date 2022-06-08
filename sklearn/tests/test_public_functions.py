from importlib import import_module
from inspect import signature

import pytest

from sklearn.utils._param_validation import generate_invalid_param_val
from sklearn.utils._param_validation import generate_valid_param
from sklearn.utils._param_validation import make_constraint


PARAM_VALIDATION_FUNCTION_LIST = [
    "sklearn.cluster.kmeans_plusplus",
]


@pytest.mark.parametrize("func_module", PARAM_VALIDATION_FUNCTION_LIST)
def test_function_param_validation(func_module):
    """Check that an informative error is raised when the value of a parameter does not
    have an appropriate type or value.
    """
    module_name, func_name = func_module.rsplit(".", 1)
    module = import_module(module_name)
    func = getattr(module, func_name)

    func_sig = signature(func)
    func_params = [
        p.name
        for p in func_sig.parameters.values()
        if p.kind not in (p.VAR_POSITIONAL, p.VAR_KEYWORD)
    ]
    parameter_constraints = getattr(func, "_skl_parameter_constraints")

    # generate valid values for the required parameters
    required_params = [
        p.name for p in func_sig.parameters.values() if p.default is p.empty
    ]
    required_params = {
        p: generate_valid_param(make_constraint(parameter_constraints[p][0]))
        for p in required_params
    }

    # check that there is a constraint for each parameter
    if func_params:
        err_msg = (
            "Mismatch between the parameter constraints defined in the "
            f"validate_params decorator and the parameters of {func_name}."
        )
        assert list(parameter_constraints.keys()) == func_params, err_msg

    # this object does not have a valid type for sure for all params
    param_with_bad_type = type("BadType", (), {})()

    for param_name in func_params:
        match = (
            rf"The '{param_name}' parameter of {func_name} must be .* Got .* instead."
        )

        # First, check that the error is raised if param doesn't match any valid type.
        with pytest.raises(ValueError, match=match):
            func(**{**required_params, param_name: param_with_bad_type})

        # Then, for constraints that are more than a type constraint, check that the
        # error is raised if param does match a valid type but does not match any valid
        # value for this type.
        constraints = parameter_constraints[param_name]
        constraints = [make_constraint(constraint) for constraint in constraints]

        for constraint in constraints:
            try:
                bad_value = generate_invalid_param_val(constraint)
            except NotImplementedError:
                continue

            with pytest.raises(ValueError, match=match):
                func(**{**required_params, param_name: bad_value})
