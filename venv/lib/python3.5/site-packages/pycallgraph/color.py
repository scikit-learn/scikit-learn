import colorsys


class ColorException(Exception):
    pass


class Color(object):

    def __init__(self, r, g, b, a=1):
        self.r = r
        self.g = g
        self.b = b
        self.a = a
        self.validate_all()

    @classmethod
    def hsv(cls, h, s, v, a=1):
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        return cls(r, g, b, a)

    def __str__(self):
        return '<Color {}>'.format(self.rgba_web())

    def validate_all(self):
        self.validate('r')
        self.validate('g')
        self.validate('b')
        self.validate('a')

    def validate(self, attr):
        v = getattr(self, attr)
        if not 0 <= v <= 1:
            raise ColorException('{} out of range 0 to 1: {}'.format(attr, v))

    @property
    def r255(self):
        return int(self.r * 255)

    @property
    def g255(self):
        return int(self.g * 255)

    @property
    def b255(self):
        return int(self.b * 255)

    @property
    def a255(self):
        return int(self.a * 255)

    def rgb_web(self):
        '''Returns a string with the RGB components as a HTML hex string.'''
        return '#{0.r255:02x}{0.g255:02x}{0.b255:02x}'.format(self)

    def rgba_web(self):
        '''Returns a string with the RGBA components as a HTML hex string.'''
        return '{0}{1.a255:02x}'.format(self.rgb_web(), self)

    def rgb_csv(self):
        '''Returns a string with the RGB components as CSV.'''
        return '{0.r255},{0.g255},{0.b255}'.format(self)
