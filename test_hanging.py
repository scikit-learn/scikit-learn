import faulthandler


def test_set_faulthandler():
    faulthandler.dump_traceback_later(10, repeat=True)
