from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import hashlib
import os

from PIL import Image as pil
from django.utils import six, lru_cache
from django.utils.six import StringIO, BytesIO

from image import settings
from image.settings import IMAGE_DEFAULT_QUALITY, IMAGE_DEFAULT_FORMAT
from image.storage import MEDIA_STORAGE, STATIC_STORAGE, IMAGE_CACHE_STORAGE

INTERNAL_CACHE_ROOT = "%s/_internal/" % settings.IMAGE_CACHE_ROOT
ALPHA_FORMATS = ["PNG"]


def power_to_rgb(value):
    if value <= 0.0031308:
        value *= 12.92
    else:
        value = 1.055 * pow(value, 0.416666666667) - 0.055
    return round(value * 255.0)


def rgb_to_power(value):
    value = float(value) / 255.0
    if value <= 0.04045:
        value /= 12.92
    else:
        value = pow((value + 0.055) / 1.055, 2.4)
    return value


def add_rgba_to_pixel(pixel, rgba, x_ammount, x_displacement):
    a = rgba[3]
    pa = pixel[3]
    if a == 1.0 and pa == 1.0:
        total_ammount = x_ammount + x_displacement
        rgba_ammount = x_ammount / total_ammount
        pixel_ammount = x_displacement / total_ammount
        return (
            pixel[0] * pixel_ammount + rgba[0] * rgba_ammount,
            pixel[1] * pixel_ammount + rgba[1] * rgba_ammount,
            pixel[2] * pixel_ammount + rgba[2] * rgba_ammount,
            pa * pixel_ammount + a * rgba_ammount,
        )
    else:
        total_ammount = x_ammount + x_displacement
        rgba_ammount = x_ammount / total_ammount
        # rgba_ammount_alpha = rgba_ammount * a
        pixel_ammount = x_displacement / total_ammount
        # pixel_ammount_alpha = pixel_ammount * pa
        return (
            pixel[0] * pixel_ammount + rgba[0] * rgba_ammount,
            pixel[1] * pixel_ammount + rgba[1] * rgba_ammount,
            pixel[2] * pixel_ammount + rgba[2] * rgba_ammount,
            pa * pixel_ammount + a * rgba_ammount,
        )


def resizeScale(img, width, height, force):
    src_width, src_height = img.size
    src_ratio = float(src_width) / float(src_height)

    if force:
        max_width = width
        max_height = height
    else:
        max_width = min(width, src_width)
        max_height = min(height, src_height)

    dst_width = max_width
    dst_height = dst_width / src_ratio

    if dst_height > max_height:
        dst_height = max_height
        dst_width = dst_height * src_ratio

    # img_width, img_height = img.size
    img = img.resize((int(dst_width), int(dst_height)), pil.ANTIALIAS)

    return img


def resizeCrop(img, width, height, center, force):
    """
    # Esto no hace nada perceptible
    RATIO = 2
    while img.size[0] / RATIO > width or img.size[1] / RATIO > height:
        img = img.resize((int(img.size[0]/RATIO), int(img.size[1]/RATIO)), pil.ANTIALIAS)
    """

    max_width = width
    max_height = height

    if not force:
        img.thumbnail((max_width, max_height), pil.ANTIALIAS)
    else:
        src_width, src_height = img.size
        src_ratio = float(src_width) / float(src_height)
        dst_width, dst_height = max_width, max_height
        dst_ratio = float(dst_width) / float(dst_height)

        if dst_ratio < src_ratio:
            crop_height = src_height
            crop_width = crop_height * dst_ratio
            x_offset = float(src_width - crop_width) / 2
            y_offset = 0
        else:
            crop_width = src_width
            crop_height = crop_width / dst_ratio
            x_offset = 0
            y_offset = float(src_height - crop_height) / 2

        center_x, center_y = center.split(',')
        center_x = float(center_x)
        center_y = float(center_y)

        x_offset = min(
            max(0, center_x * src_width - crop_width / 2),
            src_width - crop_width
        )
        y_offset = min(
            max(0, center_y * src_height - crop_height / 2),
            src_height - crop_height
        )

        img = img.crop(
            (int(x_offset), int(y_offset), int(x_offset) + int(crop_width), int(y_offset) + int(crop_height)))
        img = img.resize((int(dst_width), int(dst_height)), pil.ANTIALIAS)

    return img


def do_tint(img, tint):
    if not tint or tint is 'None':
        return

    if img.mode != "RGBA":
        img = img.convert("RGBA")

    try:
        tint_red = float(int("0x%s" % tint[0:2], 16)) / 255.0
    except ValueError:
        tint_red = 1.0

    try:
        tint_green = float(int("0x%s" % tint[2:4], 16)) / 255.0
    except ValueError:
        tint_green = 1.0

    try:
        tint_blue = float(int("0x%s" % tint[4:6], 16)) / 255.0
    except ValueError:
        tint_blue = 1.0

    try:
        tint_alpha = float(int("0x%s" % tint[6:8], 16)) / 255.0
    except ValueError:
        tint_alpha = 1.0

    try:
        intensity = float(int("0x%s" % tint[8:10], 16))
    except ValueError:
        intensity = 255.0

    if intensity > 0.0 and (tint_red != 1.0 or tint_green != 1.0 or tint_blue != 1.0 or tint_alpha != 1.0):
        # Only tint if the color provided is not ffffffff, because that equals no tint

        pixels = img.load()
        if intensity == 255.0:
            for y in range(img.size[1]):
                for x in range(img.size[0]):
                    data = pixels[x, y]
                    pixels[x, y] = (
                        int(float(data[0]) * tint_red),
                        int(float(data[1]) * tint_green),
                        int(float(data[2]) * tint_blue),
                        int(float(data[3]) * tint_alpha),
                    )
        else:
            intensity = intensity / 255.0
            intensity_inv = 1 - intensity
            tint_red *= intensity
            tint_green *= intensity
            tint_blue *= intensity
            tint_alpha *= intensity
            for y in range(img.size[1]):
                for x in range(img.size[0]):
                    data = pixels[x, y]
                    pixels[x, y] = (
                        int(float(data[0]) * intensity_inv + float(data[0]) * tint_red),
                        int(float(data[1]) * intensity_inv + float(data[1]) * tint_green),
                        int(float(data[2]) * intensity_inv + float(data[2]) * tint_blue),
                        int(float(data[3]) * intensity_inv + float(data[3]) * tint_alpha),
                    )


def do_grayscale(img):
    return img.convert('LA').convert('RGBA')


def do_paste(img, overlay, position):
    if overlay.mode != 'RGBA':
        overlay = overlay.convert('RGBA')

    # img.paste(overlay, position, overlay)
    # return img

    overlay_full = pil.new('RGBA', img.size, color=(255, 255, 255, 0))
    overlay_full.paste(overlay, position)

    # img.paste(overlay_full, (0,0), overlay_full)
    # return img

    return pil.alpha_composite(img, overlay_full)

    # overlay_pixels = overlay.load()
    # img_pixels = img.load()
    # overlay_width, overlay_height = overlay.size
    # x_offset, y_offset = position
    #
    # for y in range(min(overlay_height, img.size[1] - y_offset)):
    #     for x in range(min(overlay_width, img.size[0] - x_offset)):
    #         img_pixel = img_pixels[x + x_offset, y + y_offset]
    #         overlay_pixel = overlay_pixels[x, y]
    #         ia = img_pixel[3]
    #         oa = overlay_pixel[3]
    #         if oa == 0:
    #             # overlay is transparent, nothing to do
    #             continue
    #         elif oa == 255:
    #             # overlay is opaque, ignore img pixel
    #             new_pixel = overlay_pixel
    #         elif ia == 0:
    #             # image pixel is 100% transparent, only overlay matters
    #             new_pixel = overlay_pixel
    #         elif ia == 255:
    #             # simpler math
    #             oa = float(oa) / 255.0
    #             oa1 = 1.0 - oa
    #             new_pixel = (
    #                          int(power_to_rgb( rgb_to_power(img_pixel[0]) * oa1 + rgb_to_power(overlay_pixel[0]) * oa  )),
    #                          int(power_to_rgb( rgb_to_power(img_pixel[1]) * oa1 + rgb_to_power(overlay_pixel[1]) * oa  )),
    #                          int(power_to_rgb( rgb_to_power(img_pixel[2]) * oa1 + rgb_to_power(overlay_pixel[2]) * oa  )),
    #                          255,
    #                          )
    #         else:
    #             # complex math
    #             oa = float(oa) / 255.0
    #             ia = float(ia) / 255.0
    #             oa1 = 1 - oa
    #             #total_alpha_percent = oa + ia
    #             overlay_percent = oa
    #             image_percent = ia * oa1
    #
    #             new_pixel = (
    #                          int(power_to_rgb( rgb_to_power(img_pixel[0]) * image_percent + rgb_to_power(overlay_pixel[0]) * overlay_percent  )),
    #                          int(power_to_rgb( rgb_to_power(img_pixel[1]) * image_percent + rgb_to_power(overlay_pixel[1]) * overlay_percent  )),
    #                          int(power_to_rgb( rgb_to_power(img_pixel[2]) * image_percent + rgb_to_power(overlay_pixel[2]) * overlay_percent  )),
    #                          int((oa + ia * oa1) * 255.0),
    #                          )
    #
    #         img_pixels[x + x_offset, y + y_offset] = new_pixel


def do_overlay(img, overlay_path, overlay_source=None, overlay_tint=None, overlay_size=None, overlay_position=None):
    if not overlay_path:
        return img

    if overlay_source == 'media':
        overlay = pil.open(MEDIA_STORAGE.open(overlay_path))
    else:
        overlay = pil.open(STATIC_STORAGE.open(overlay_path))

    # We want the overlay to fit in the image
    iw, ih = img.size
    ow, oh = overlay.size
    overlay_ratio = float(ow) / float(oh)

    if overlay_size:
        tw, th = overlay_size.split(',')
        ow = int(round(float(tw.strip()) * iw))
        oh = int(round(float(th.strip()) * ih))
        if ow < 0:
            ow = oh * overlay_ratio
        elif oh < 0:
            oh = ow / overlay_ratio

        overlay = resizeScale(overlay, ow, oh, overlay_source + "/" + overlay_path)
        ow, oh = overlay.size
    else:
        have_to_scale = False
        if ow > iw:
            ow = iw
            oh = int(float(iw) / overlay_ratio)
            have_to_scale = True
        if oh > ih:
            ow = int(float(ih) * overlay_ratio)
            oh = ih
            have_to_scale = True

        if have_to_scale:
            overlay = resizeScale(overlay, ow, oh, overlay_source + "/" + overlay_path)
            ow, oh = overlay.size

    if overlay_tint:
        do_tint(overlay, overlay_tint)

    if not overlay_position:
        target_x = int((iw - ow) / 2)
        target_y = int((ih - oh) / 2)
    else:
        tx, ty = overlay_position.split(',')
        tx = tx.strip()
        ty = ty.strip()

        if tx == "":
            # Center horizontally.
            target_x = int((iw - ow) / 2)
        else:
            if "!" in tx:
                # X origin on the right side.
                x_percent = 1.0 - float(tx.replace("!", "").strip())
                target_x = int(round(x_percent * iw) - ow)
            else:
                # X origin on the left side.
                x_percent = float(tx)
                target_x = int(round(x_percent * iw))

        if ty == "":
            # Center vertically.
            target_y = int((ih - oh) / 2)
        else:
            if "!" in ty:
                # Y origin on the bottom side.
                y_percent = 1.0 - float(ty.replace("!", "").strip())
                target_y = int(round(y_percent * ih) - oh)
            else:
                # Y origin on the top side.
                y_percent = float(ty)
                target_y = int(round(y_percent * ih))

    """
    TODO: paste seems to be buggy, because pasting over opaque background returns a non opaque image
    (the parts that are not 100% opaque or 100% transparent become partially transparent. 
    the putalpha workareound doesn't seem to look nice enough
    """
    img = do_paste(img, overlay, (target_x, target_y))

    return img


def do_overlays(img, overlays, overlay_tints, overlay_sources, overlay_sizes, overlay_positions):
    overlay_index = 0

    for overlay in overlays:

        try:
            overlay_tint = overlay_tints[overlay_index]
        except (IndexError, TypeError):
            overlay_tint = None

        if overlay_tint == "None":
            overlay_tint = None

        try:
            overlay_source = overlay_sources[overlay_index]
        except (IndexError, TypeError):
            overlay_source = 'static'

        try:
            overlay_size = overlay_sizes[overlay_index]
        except (IndexError, TypeError):
            overlay_size = None

        if overlay_size == "None":
            overlay_size = None

        try:
            overlay_position = overlay_positions[overlay_index]
        except (IndexError, TypeError):
            overlay_position = None

        if overlay_position == "None":
            overlay_position = None

        img = do_overlay(img, overlay, overlay_source, overlay_tint, overlay_size, overlay_position)
        overlay_index += 1

    return img


def do_mask(img, mask_path, mask_source, mask_mode=None):
    if not mask_path:
        return img

    if mask_source == 'media':
        mask = pil.open(MEDIA_STORAGE.open(mask_path)).convert("RGBA")
    else:
        mask = pil.open(STATIC_STORAGE.open(mask_path)).convert("RGBA")

    # We want the mask to have the same size than the image
    if mask_mode == 'distort':
        iw, ih = img.size
        mw, mh = mask.size
        if mw != iw or mh != ih:
            mask = mask.resize((iw, ih), pil.ANTIALIAS)

    else:
        # We want the overlay to fit in the image
        iw, ih = img.size
        ow, oh = mask.size
        overlay_ratio = float(ow) / float(oh)
        have_to_scale = False
        if ow > iw:
            ow = iw
            oh = int(float(iw) / overlay_ratio)
        if oh > ih:
            ow = int(float(ih) * overlay_ratio)
            oh = ih

        if ow != iw or oh != ih:
            have_to_scale = True

        if have_to_scale:
            nmask = mask.resize((ow, oh), pil.ANTIALIAS)
            mask = pil.new('RGBA', (iw, ih))
            # mask.paste(nmask, (int((iw - ow) / 2), int((ih - oh) / 2)), nmask)
            mask = do_paste(mask, nmask, (int((iw - ow) / 2), int((ih - oh) / 2)))
            ow, oh = mask.size

    r, g, b, a = mask.split()
    img.putalpha(a)


def do_fill(img, fill, width, height):
    if not fill:
        return img

    overlay = img

    fill_color = (
        int("0x%s" % fill[0:2], 16),
        int("0x%s" % fill[2:4], 16),
        int("0x%s" % fill[4:6], 16),
        int("0x%s" % fill[6:8], 16),
    )
    img = pil.new("RGBA", (width, height), fill_color)

    iw, ih = img.size
    ow, oh = overlay.size

    # img.paste(overlay, (int((iw - ow) / 2), int((ih - oh) / 2)), overlay)
    img = do_paste(img, overlay, (int((iw - ow) / 2), int((ih - oh) / 2)))

    return img


def do_padding(img, padding):
    if not padding:
        return img
    try:
        padding = float(padding) * 2.0
        if padding > .9:
            padding = .9
        if padding <= 0.0:
            return img
    except ValueError:
        return

    iw, ih = img.size

    img.thumbnail(
        (
            int(round(float(img.size[0]) * (1.0 - padding))),
            int(round(float(img.size[1]) * (1.0 - padding)))
        ),
        pil.ANTIALIAS
    )

    img = do_fill(img, "ffffff00", iw, ih)

    return img


def do_background(img, background):
    if not background:
        return img

    overlay = img

    fill_color = (
        int("0x%s" % background[0:2], 16),
        int("0x%s" % background[2:4], 16),
        int("0x%s" % background[4:6], 16),
        int("0x%s" % background[6:8], 16),
    )
    img = pil.new("RGBA", overlay.size, fill_color)

    iw, ih = img.size
    ow, oh = overlay.size

    # img.paste(overlay, (int((iw - ow) / 2), int((ih - oh) / 2)), overlay)
    img = do_paste(img, overlay, (int((iw - ow) / 2), int((ih - oh) / 2)))

    return img


def do_rotate(img, rotation):
    if not rotation:
        return img

    try:
        rotation = float(rotation)
    except ValueError:
        rotation = 0

    if rotation % 90 == 0:
        img = img.rotate(rotation, pil.NEAREST, expand=True)
    else:
        img = img.rotate(rotation, pil.BICUBIC, expand=True)

    return img


def render(data, width, height, force=True, padding=None, overlays=(), overlay_sources=(),
           overlay_tints=(), overlay_sizes=None, overlay_positions=None, mask=None, mask_source=None,
           center=".5,.5", format=IMAGE_DEFAULT_FORMAT, quality=IMAGE_DEFAULT_QUALITY, fill=None, background=None,
           tint=None, pre_rotation=None, post_rotation=None, crop=True, grayscale=False):
    """
    Rescale the given image, optionally cropping it to make sure the result image has the specified width and height.
    """

    if not isinstance(data, six.string_types):
        input_file = BytesIO(data)
    else:
        input_file = StringIO(data)

    img = pil.open(input_file)
    if img.mode != "RGBA":
        img = img.convert("RGBA")

    if width is None:
        width = img.size[0]
    if height is None:
        height = img.size[1]

    img = do_rotate(img, pre_rotation)

    if crop:
        img = resizeCrop(img, width, height, center, force)
    else:
        img = resizeScale(img, width, height, force)

    if grayscale:
        img = do_grayscale(img)
    do_tint(img, tint)
    img = do_fill(img, fill, width, height)
    img = do_background(img, background)
    do_mask(img, mask, mask_source)
    img = do_overlays(img, overlays, overlay_tints, overlay_sources, overlay_sizes, overlay_positions)
    img = do_padding(img, padding)
    img = do_rotate(img, post_rotation)

    tmp = BytesIO()

    if not format.upper() in ALPHA_FORMATS:
        img = img.convert("RGB")

    img.save(tmp, format, quality=quality)
    tmp.seek(0)
    output_data = tmp.getvalue()
    input_file.close()
    tmp.close()

    return output_data


@lru_cache.lru_cache(maxsize=128)
def image_create_token(parameters):
    return "image_token_%s" % hashlib.sha1(parameters.encode("utf8")).hexdigest()


def image_tokenize(session, parameters):
    if session:
        token = None
        for k, v in session.items():
            if v == parameters:
                token = k
                break
        if token is None:
            token = image_create_token(parameters)
            session[token] = parameters
    else:
        token = image_create_token(parameters)
    return token


def image_url(session, parameters, image_field, generate=False):
    if generate:
        from image import views as image_views

        autogen = 'autogen=true' in parameters
        image_views.image(session, str(image_field), parameters, autogen)

    image_path = os.path.join(image_tokenize(session, parameters), six.text_type(image_field))
    return IMAGE_CACHE_STORAGE.url(image_path)
