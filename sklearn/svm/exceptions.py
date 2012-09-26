
class SolverTimeout(Exception):
    def __init__(self, fit_rval):
        self.fit_rval = fit_rval


