# -*- coding: UTF-8 -*-

from django import template
from django.db.models.fields.files import ImageFieldFile
from django.http.request import HttpRequest
from django.template.defaulttags import register
from django.utils import six

from image import views as image_views
from image.storage import IMAGE_CACHE_STORAGE
from image.utils import image_url
from image.video_field import VideoFieldFile


class ImageNode(template.Node):
    def __init__(self, image_field, parameters):
        self.image_field = image_field
        self.parameters = parameters  # It is not a variable because we compiled filters.

    def render(self, context):
        try:
            request = context['request']
            session = request.session
        except KeyError:
            request = HttpRequest()
            session = None

        image_field = self.image_field.resolve(context)
        try:
            parameters = self.parameters.resolve(context)
        except template.VariableDoesNotExist:
            parameters = self.parameters

        if isinstance(image_field, VideoFieldFile):
            parameters += "&video=true"

        if isinstance(image_field, ImageFieldFile) or isinstance(image_field, VideoFieldFile):
            try:
                parameters = parameters + "&center=" + image_field.__image_center_instance__.__unicode__()
            except AttributeError:
                pass

        if "autogen=true" in parameters or getattr(IMAGE_CACHE_STORAGE, "autogen_required", False):
            # We want the image to be generated immediately
            image_views.image(request, six.text_type(image_field), parameters, True)

        # image_path = os.path.join(image_tokenize(session, parameters), six.text_type(image_field))
        # return IMAGE_CACHE_STORAGE.url(image_path)
        return image_url(session, parameters, image_field)

        # return reverse(
        #     'image.views.image',
        #     args=(
        #         str(image_field),
        #         image_tokenize(session, parameters)
        #     )
        # )


@register.tag(name="image")
def image(parser, token):
    try:
        # split_contents() knows not to split quoted strings.
        tag_name, image_field, parameters = token.split_contents()
    except ValueError:
        raise template.TemplateSyntaxError("%r tag requires 2 arguments " % token.contents.split()[0])

    image_field = parser.compile_filter(image_field)
    parameters = parser.compile_filter(parameters)

    return ImageNode(image_field, parameters)
