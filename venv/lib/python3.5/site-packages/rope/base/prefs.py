class Prefs(object):

    def __init__(self):
        self.prefs = {}
        self.callbacks = {}

    def set(self, key, value):
        """Set the value of `key` preference to `value`."""
        if key in self.callbacks:
            self.callbacks[key](value)
        else:
            self.prefs[key] = value

    def add(self, key, value):
        """Add an entry to a list preference

        Add `value` to the list of entries for the `key` preference.

        """
        if not key in self.prefs:
            self.prefs[key] = []
        self.prefs[key].append(value)

    def get(self, key, default=None):
        """Get the value of the key preference"""
        return self.prefs.get(key, default)

    def add_callback(self, key, callback):
        """Add `key` preference with `callback` function

        Whenever `key` is set the callback is called with the
        given `value` as parameter.

        """
        self.callbacks[key] = callback

    def __setitem__(self, key, value):
        self.set(key, value)

    def __getitem__(self, key):
        return self.get(key)
