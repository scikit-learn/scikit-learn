class _JediError(Exception):
    pass


class InternalError(_JediError):
    pass


class WrongVersion(_JediError):
    pass
