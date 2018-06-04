# -*- coding: UTF-8 -*-
from django.db import models
from django.db.models.signals import post_init
from django.db.models.fields import FieldDoesNotExist
from django.utils import six

from image.video_field import VideoField
from image.forms import ImageCenterFormField


class ImageCenter(object):
    def __init__(self, image_field, x=None, y=None, xy=None):
        self.image_field = image_field
        if (not x is None and y is None) or (x is None and not y is None):
            raise ValueError(u"If x or y are provided, both have to be provided: x=" + six.text_type(x) + u", y=" + six.text_type(y))
        if not x is None and not y is None:
            if x < 0 or x > 1 or y < 0 or y > 1:
                raise ValueError(u"Valid values for x and y go from 0 to 1: x=" + six.text_type(x) + u", y=" + six.text_type(y))
            self.x = float(x)
            self.y = float(y)
        else:
            if xy is None:
                self.x = .5
                self.y = .5
            else:
                try:
                    x, y = xy.split(",")
                except ValueError:
                    x = .5
                    y = .5
                self.x = float(x)
                self.y = float(y)

    def __unicode__(self):
        return str(self.x) + "," + str(self.y)


class ImageCenterField(models.Field):

    attr_class = ImageCenter

    description = "A field that stores the center of attention for an image."

    # __metaclass__ = models.SubfieldBase

    def __init__(self, image_field=None, *args, **kwargs):
        if image_field is not None:
            if not isinstance(image_field, models.ImageField) and  not isinstance(image_field, VideoField):
                raise ValueError("image_field value must be an ImageField or VideoField instance")
        kwargs["default"] = ".5,.5"
        self.image_field = image_field
        super(ImageCenterField, self).__init__(*args, **kwargs)

    def set_instance(self, instance):
        self.instance = instance

    def formfield(self, **kwargs):
        defaults = {'form_class': ImageCenterFormField}
        defaults.update(kwargs)
        return super(ImageCenterField, self).formfield(**defaults)

    def db_type(self, connection):
        return "char(100)"

    # Esta función es llamada al leer un valor de la base de datos
    def to_python(self, value):
        if isinstance(value, ImageCenter):
            return value

        return ImageCenter(self.image_field, xy=value)

    # Esta función es llamada al escribir un valor en la base de datos
    def get_db_prep_value(self, value, connection=None, prepared=False):
        try:
            return str(value.x) + "," + str(value.y)
        except AttributeError:
            return str(value)

    def from_db_value(self, value, expression, connection, context):
        return ImageCenter(self.image_field, xy=value)

    def value_to_string(self, obj):
        value = getattr(obj, self.attname)
        return self.get_db_prep_value(value)

    def query_string(self):
        return "center=" + str(self.x) + "," + str(self.y)


def post_init_capture(sender, instance, *args, **kwargs):
    fields = instance.__class__._meta.get_fields()
    for field in fields:
        if isinstance(field, ImageCenterField):
            image_field = instance.__class__._meta.get_field(field.image_field.name)
            image_instance = instance.__getattribute__(image_field.name)
            image_center_instance = instance.__getattribute__(field.name)
            image_instance.__image_center_instance__ = image_center_instance
            if isinstance(image_center_instance, six.string_types):
                image_center_instance = ImageCenter(image_field, xy=image_center_instance)
                setattr(instance, field.name, image_center_instance)
            image_center_instance.image_path = six.text_type(image_instance)

post_init.connect(post_init_capture)

try:
    from south.modelsinspector import add_introspection_rules
    add_introspection_rules([], ["^image\.fields\.ImageCenterField$"])
    add_introspection_rules([], ["^image\.video_field\.VideoField$"])
except ImportError:
    pass
