# -*- coding: UTF-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from django import forms
from django.utils.safestring import mark_safe
from django.utils import six
from django.forms.utils import flatatt
from django.utils.encoding import force_text
try:
    from django.core.urlresolvers import reverse
except ImportError:
    from django.urls import reverse
import threading


from image.video_field import VideoField


COUNTER = 0


class ImageCenterFormWidget(forms.Widget):

    #def __init__(self, attrs=None):
    #    super(ImageCenterFormWidget, self).__init__(attrs)

    def _format_value(self, value):
        return six.text_type(value)

    def render(self, name, value, attrs=None, **kwargs):

        global COUNTER

        if value is None:
            value = ''
        try:
            final_attrs = self.build_attrs(attrs, name=name)
        except TypeError:
            final_attrs = self.build_attrs(attrs, {"name": name})
        if value != '':
            # Only add the 'value' attribute if a value is non-empty.
            final_attrs['value'] = force_text(self._format_value(value))

        resp = ''
        if getattr(value, 'image_path', None):
            try:
                extra_parms = ""
                if isinstance(value.image_field, VideoField):
                    extra_parms += "&video=true"

                extra_parms += "&is_admin=true"

                resp = '<div style="display:inline-block; position:relative; border:1px solid black;">'
                resp += '<img id="image_center-' + str(COUNTER) + '" src="' + reverse(
                    'image.views.image', args=('format=png&width=150&height=150&mode=scale' + extra_parms, value.image_path)) + '" onclick=""/>'
                resp += '<img id="image_center_crosshair-' + str(COUNTER) + '" src="' + reverse('image.views.crosshair') + '" style="position:absolute; left:0; top:0;" />'
                resp += '</div>'
                resp += '<script>'
                resp += '(function($) {'
                resp += '    $(window).load(function(){'
                resp += '        var crosshair = document.getElementById("image_center_crosshair-' + str(COUNTER) + '");'
                resp += '        var image = document.getElementById("image_center-' + str(COUNTER) + '");'
                resp += '        var iw = $(image).width();'
                resp += '        var ih = $(image).height();'
                resp += '        $(crosshair).css( { left : (iw*' + str(value.x) + ' - 7)+"px", top : (ih*' + str(value.y) + ' - 7)+"px" } );'
                resp += '        $(image).parent().parent().find("input").hide();'
                resp += '        $(image).parent().click(function(e){'
                resp += '            var nx = e.pageX - $(image).offset().left;'
                resp += '            var ny = e.pageY - $(image).offset().top;'
                resp += '            crosshair.style.left=(nx - 7)+"px";'
                resp += '            crosshair.style.top=(ny - 7)+"px";'
                resp += '            $(image).parent().parent().find("input").val( (nx/iw)+","+(ny/ih) );'
                resp += '        });'
                resp += '});'
                resp += '})(django.jQuery);'
                resp += '</script>'
                resp += u'<input%s />' % flatatt(final_attrs)

                lock = threading.Lock()
                with lock:
                    COUNTER += 1
                    if COUNTER > 4000000000:
                        COUNTER = 0
            except AttributeError:
                resp = 'Only available once the image has been saved'

        return mark_safe(resp)


class ImageCenterFormField(forms.Field):

    widget = ImageCenterFormWidget

    def __init__(self, *args, **kwargs):
        kwargs['required'] = False
        super(ImageCenterFormField, self).__init__(*args, **kwargs)

    def clean(self, value):
        value = self.to_python(value)
        return value
