from django.db import models
from django.db.models.fields.files import FieldFile

from django.core.files import File


def get_video_dimensions(path):
    from ffvideo import VideoStream

    vs = VideoStream(path)
    return (vs.frame_width, vs.frame_height)


class VideoFile(File):
    """
    A mixin for use alongside django.core.files.base.File, which provides
    additional features for dealing with images.
    """

    def _get_width(self):
        return self._get_video_dimensions()[0]

    width = property(_get_width)

    def _get_height(self):
        return self._get_video_dimensions()[1]

    height = property(_get_height)

    def _get_video_dimensions(self):
        if not hasattr(self, '_dimensions_cache'):
            close = self.closed
            self.open()
            self._dimensions_cache = get_video_dimensions(self.path)
        return self._dimensions_cache


# A video field is exactly a file field with a different signature
class VideoFieldFile(VideoFile, FieldFile):
    pass


class VideoField(models.FileField):
    attr_class = VideoFieldFile
