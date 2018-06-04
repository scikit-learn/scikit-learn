#-*- coding: UTF-8 -*-

try:
    from django.core.urlresolvers import reverse
except ImportError:
    from django.urls import reverse
from django.utils import six

from image.views import image as image_view
from image.utils import image_create_token


def get_image_url(image, parameters):
    if 'autogen=true' in parameters:
        image_view(None, str(image), parameters, True)
    
    return reverse(
        'image.views.image',
        args=(
            image_create_token(parameters),
            six.text_type(image),
        )
    )


