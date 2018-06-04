# -*- coding: UTF-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import hashlib

import requests
import six
from django.core.files.base import ContentFile
from encodings.base64_codec import base64_decode
import os
import urllib
import traceback

from django.http import HttpResponse, QueryDict
from django.http.response import Http404
from django.utils import timezone
from django.utils.encoding import force_text

from image.settings import IMAGE_CACHE_HTTP_EXPIRATION, IMAGE_CACHE_ROOT
from image.storage import IMAGE_CACHE_STORAGE, MEDIA_STORAGE, STATIC_STORAGE
from image.utils import IMAGE_DEFAULT_FORMAT, IMAGE_DEFAULT_QUALITY,\
    image_create_token, render
from image.videothumbs import generate_thumb


def image(request, path, token, autogen=False):

    original_token = token
    token = original_token.split("-:-")[0]

    has_admin_perm = request.user.has_perm('admin') if request else False
    is_admin = False
    if ("is_admin=true" in token and request and has_admin_perm) or autogen:
        parameters = token
        is_admin = True
        if autogen:
            token = image_create_token(parameters)
    else:
        parameters = request.session.get(token, token)

    cached_image_file = os.path.join(path, token)

    now = timezone.now()
    expire_offset = timezone.timedelta(seconds=IMAGE_CACHE_HTTP_EXPIRATION)

    response = HttpResponse()
    response['Content-type'] = 'image/jpeg'
    response['Expires'] = (now + expire_offset).strftime("%a, %d %b %Y %T GMT")
    response['Last-Modified'] = now.strftime("%a, %d %b %Y %T GMT")
    response['Cache-Control'] = 'max-age=3600, must-revalidate'
    response.status_code = 200

    # If we already have the cache we send it instead of recreating it
    if IMAGE_CACHE_STORAGE.exists(cached_image_file):
        
        if autogen:
            return 'Already generated'
        
        try:
            f = IMAGE_CACHE_STORAGE.open(cached_image_file, "rb")
        except IOError:
            raise Http404()
        response.write(f.read())
        f.close()

        try:
            # Django 2.0 support
            modified_time = IMAGE_CACHE_STORAGE.get_modified_time(cached_image_file)
        except AttributeError:
            modified_time = IMAGE_CACHE_STORAGE.modified_time(cached_image_file)
        response['Last-Modified'] = modified_time.strftime("%a, %d %b %Y %T GMT")
        return response
    
    if parameters == token and not is_admin:
        if has_admin_perm and 'is_admin=true' in path:
            #  Retrocompatibilidad.
            return image(request, original_token, path, autogen=autogen)
        return HttpResponse("Forbidden", status=403)

    qs = QueryDict(parameters)

    file_storage = MEDIA_STORAGE
    if qs.get('static', '') == "true":
        file_storage = STATIC_STORAGE

    format = qs.get('format', IMAGE_DEFAULT_FORMAT)
    quality = int(qs.get('quality', IMAGE_DEFAULT_QUALITY))
    mask = qs.get('mask', None)
    mask_source = qs.get('mask_source', None)

    if mask is not None:
        format = "PNG"

    fill = qs.get('fill', None)
    background = qs.get('background', None)
    tint = qs.get('tint', None)

    center = qs.get('center', ".5,.5")
    mode = qs.get('mode', "crop")
    enlarge = qs.get('enlarge', "true")

    overlays = qs.getlist('overlay')
    overlay_sources = qs.getlist('overlay_source')
    overlay_tints = qs.getlist('overlay_tint')
    overlay_sizes = qs.getlist('overlay_size')
    overlay_positions = qs.getlist('overlay_position')

    width = qs.get('width', None)
    if width:
        width = int(width)
    height = qs.get('height', None)
    if height:
        height = int(height)

    pre_rotation = qs.get('pre_rotation', None)
    post_rotation = qs.get('post_rotation', None)

    try:
        padding = float(qs.get('padding', None))
    except TypeError:
        padding = 0.0

    grayscale = bool(qs.get('grayscale', False))

    if "video" in qs:
        data, http_response = generate_thumb(file_storage, force_text(path))
        response.status_code = http_response
    else:
        try:
            try:
                f = requests.get(qs['url'])
                data = f.content
                f.close()
            except KeyError:
                f = file_storage.open(path)
                data = f.read()
                f.close()
        except IOError:
            response.status_code = 404
            data = ""

    if data:
        try:
            crop = (mode != "scale")
            force = (enlarge == "true")
            output_data = render(data, width, height, force=force, padding=padding, overlays=overlays,
                                 overlay_sources=overlay_sources, overlay_tints=overlay_tints,
                                 overlay_positions=overlay_positions, overlay_sizes=overlay_sizes, mask=mask,
                                 mask_source=mask_source, center=center, format=format, quality=quality, fill=fill,
                                 background=background, tint=tint, pre_rotation=pre_rotation,
                                 post_rotation=post_rotation, crop=crop, grayscale=grayscale)
        except IOError:
            traceback.print_exc()
            response.status_code = 500
            output_data = ""
    else:
        output_data = data

    if response.status_code == 200:
        IMAGE_CACHE_STORAGE.save(cached_image_file,  output_data, )
        if autogen:
            return 'Generated ' + six.text_type(response.status_code)
    else:
        if autogen:
            return 'Failed ' + cached_image_file
    
    response.write(output_data)

    return response


def crosshair(request):

    response = HttpResponse()
    response['Content-type'] = 'image/png'
    response['Expires'] = 'Fri, 09 Dec 2327 08:34:31 GMT'
    response['Last-Modified'] = 'Fri, 24 Sep 2010 11:36:29 GMT'
    output, length = base64_decode(b'iVBORw0KGgoAAAANSUhEUgAAAA8AAAAPCAYAAAA71pVKAAAAwElEQVQoz6WTwY0CMQxFX7DEVjCShUQnc6YCdzCH1EYDboICphb28veA2UULSBHzLpEif9vfcRr/kHQF9jzz3Vr74hWSLpKUmYoIubvMTO6uiFBmqri8FPbeBazAAhwBq3MB1t77c4IH4flNy9T9+Z7g12Nm3iu+Ez4mWMvCFUmKCFVrIywRcasuSe6u8jbC0d3/xGamGs4IZmaSpB3ANE0Ah0HxoeLZAczzDHAaFJ8qfuO0N73z5g37dLfbll/1A+4O0Wm4+ZiPAAAAAElFTkSuQmCC')
    response.write(output)
    return response
