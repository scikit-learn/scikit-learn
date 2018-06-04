from django.core.files.storage import FileSystemStorage
import os
from django.db.models.signals import post_save, post_delete, pre_save
from django.db.models.fields.files import FileField
from django.core.exceptions import ObjectDoesNotExist, MultipleObjectsReturned

from django.conf import settings as django_settings
from image.storage import IMAGE_CACHE_STORAGE


def safe_delete(path):
    if isinstance(IMAGE_CACHE_STORAGE, FileSystemStorage):
        full_path = os.path.join(IMAGE_CACHE_STORAGE.location, path)
        if os.path.isdir(full_path):
            os.rmdir(full_path)
            return
    IMAGE_CACHE_STORAGE.delete(path)


def remove_directory(dir_path):
    try:
        # Since no all storages support exists for directories, we check for OSError
        contents = IMAGE_CACHE_STORAGE.listdir(dir_path)
    except OSError:
        pass
    else:
        for directory in contents[0]:
            safe_delete(os.path.join(dir_path, directory))
        for filename in contents[1]:
            safe_delete(os.path.join(dir_path, filename))

    if IMAGE_CACHE_STORAGE.exists(dir_path):
        # In some storages like amazon S3 there are no directories
        safe_delete(dir_path)


def remove_cache(image_path):
    if image_path:
        remove_directory(image_path)


def prepare_image_cache_cleanup(sender, instance=None, **kwargs):
    if instance is None:
        return
    instance.old_image_fields = {}

    old_instance = None

    for field in instance._meta.fields:
        if isinstance(field, FileField):
            if not old_instance:
                try:
                    old_instance = sender.objects.get(pk=instance.pk)
                except (ObjectDoesNotExist, MultipleObjectsReturned):
                    return

            instance.old_image_fields[field.attname] = field.value_to_string(old_instance)


def clear_prepared_image_cache_cleanup(sender, instance=None, created=False, **kwargs):
    if created:
        return
    if instance is None:
        return
    for field in instance._meta.fields:
        if isinstance(field, FileField):
            if instance.old_image_fields[field.attname] != field.value_to_string(instance):
                remove_cache(instance.old_image_fields[field.attname])


def clear_image_cache(sender, instance, **kwargs):
    for field in instance._meta.fields:
        if isinstance(field, FileField):
            remove_cache(field.value_to_string(instance))


pre_save.connect(prepare_image_cache_cleanup)
post_save.connect(clear_prepared_image_cache_cleanup)
post_delete.connect(clear_image_cache)

#reversion compatibility
if 'reversion' in django_settings.INSTALLED_APPS:
    try:
        from reversion.models import pre_revision_commit, post_revision_commit

        pre_revision_commit.connect(prepare_image_cache_cleanup)
        post_revision_commit.connect(clear_prepared_image_cache_cleanup)
    except ImportError:
        pass

# http://bakery-app01-static.s3.amazonaws.com/image/actuality/file0001116000079.jpg/image_token_73ff82d8ce1577a8a22f5a7d29a772a3ffc6e76c

remove_directory('actuality/file0001116000079.jpg')