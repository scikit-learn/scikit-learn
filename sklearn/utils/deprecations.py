import warnings


def _deprecate_sample_weight(sample_weight, sample_props):
    if sample_props is None:
        sample_props = {}
    if "sample_weight" in sample_props.keys() and sample_weight is not None:
        raise ValueError("sample_weight and sample_props['sample_weight'] passed to fit. "
                         "Please specify only one of the two.")
    if sample_weight is not None:
        warnings.warn("The sample_weight parameter was removed by the "
                      "sample_props parameter and will be removed in 0.19.",
                      DeprecationWarning)
        return sample_weight
    else:
        return sample_props.get('sample_weight', None)
