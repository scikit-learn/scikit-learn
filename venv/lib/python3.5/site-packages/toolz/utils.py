def raises(err, lamda):
    try:
        lamda()
        return False
    except err:
        return True


no_default = '__no__default__'
