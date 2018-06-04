from jedi.evaluate.cache import evaluator_function_cache


@evaluator_function_cache()
def get_yield_exprs(evaluator, funcdef):
    return list(funcdef.iter_yield_exprs())
