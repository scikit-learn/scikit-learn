# coding=UTF-8
from django.conf import global_settings

__author__ = 'franki'

from django.conf import settings


def get(key, default):
    return getattr(settings, key, default)

IMAGE_DEFAULT_FORMAT = get('IMAGE_DEFAULT_FORMAT', 'JPEG')
IMAGE_DEFAULT_QUALITY = get('IMAGE_DEFAULT_QUALITY', 85)

IMAGE_CACHE_STORAGE = get('IMAGE_CACHE_STORAGE', 'image.storage.ImageCacheStorage')
IMAGE_CACHE_ROOT = get('IMAGE_CACHE_ROOT', None)

# If IMAGE_CACHE_URL differs from the url of the image view then you must always use autogen or have proper rewrite rules on your server
IMAGE_CACHE_URL = get('IMAGE_CACHE_URL', '/image/')


IMAGE_CACHE_HTTP_EXPIRATION = get('IMAGE_CACHE_HTTP_EXPIRATION', 3600 * 24 * 30)

FILE_UPLOAD_TEMP_DIR = get('FILE_UPLOAD_TEMP_DIR', global_settings.FILE_UPLOAD_TEMP_DIR)
STATICFILES_STORAGE = get('STATICFILES_STORAGE', global_settings.STATICFILES_STORAGE)
