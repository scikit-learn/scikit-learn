# -*- coding: UTF-8 -*-
from django.conf.urls import url

from image import views

urlpatterns = [
    url(r'^crosshair$', views.crosshair, name="image.views.crosshair"),
    url(r'^(?P<token>[:\-\w_=&]+)/(?P<path>.+)$', views.image, name="image.views.image"),
]
