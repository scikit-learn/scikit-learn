# coding=utf-8
from django.conf import settings
import sys

__author__ = 'franki'

from django.apps import AppConfig


class ImageConfig(AppConfig):
    name = 'image'
    label = 'image'
    verbose_name = "Image"

    # def ready(self):
    #     if settings.DEBUG and not 'django.template.context_processors.request' in settings.TEMPLATE_CONTEXT_PROCESSORS:
    #         print >> sys.stderr, \
    #             "image: Add 'django.template.context_processors.request' to TEMPLATE_CONTEXT_PROCESSORS in order to\n" \
    #             "give access to sessions from templates. Otherwise set autogen=true in all uses. This message only\n" \
    #             "appears with DEBUG enabled."
