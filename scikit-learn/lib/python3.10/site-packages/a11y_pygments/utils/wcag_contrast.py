# Methods to calculate WCAG contrast ratio and check if it passes AA or AAA
# https://www.w3.org/WAI/WCAG21/Understanding/contrast-minimum.html

import re

from typing import NewType, Tuple, TypeAlias, Union


# float01 is a float greater than or equal to 0 and less than or equal to 1
float01 = NewType("float01", float)
RGBColor: TypeAlias = Tuple[float01, float01, float01]


def hexdd_to_float01(xx: str) -> float01:
    """Convert a two-digit, 8-bit hex value to a float between 0 and 1.

    >>> hexdd_to_float01("00")
    0.0
    >>> hexdd_to_float01("33")
    0.2
    >>> hexdd_to_float01("55")
    0.3333333333333333
    >>> hexdd_to_float01("cc")
    0.8
    >>> hexdd_to_float01("ff")
    1.0
    """
    return float01(int(xx, 16) / 255)


def hex_to_rgb01(hex: str) -> RGBColor:
    """Convert a hex defined colour to RGB.

    Args:
        hex (string): color in hex format (#rrggbb/#rgb)

    Returns:
        rgb: tuple of rgb values ``(r, g, b)``, where each channel (red, green,
        blue) can assume values between 0 and 1 (inclusive).
    """

    if re.match(r"\A#[a-fA-F0-9]{6}\Z", hex):
        # Full hex color (#rrggbb) format
        return (
            hexdd_to_float01(hex[1:3]),
            hexdd_to_float01(hex[3:5]),
            hexdd_to_float01(hex[5:7]),
        )

    if re.match(r"\A#[a-fA-F0-9]{3}\Z", hex):
        # Short hex color (#rgb) format, shorthand for #rrggbb
        return (
            hexdd_to_float01(hex[1] * 2),
            hexdd_to_float01(hex[2] * 2),
            hexdd_to_float01(hex[3] * 2),
        )

    raise ValueError("Invalid hex color format")


def sRGB_channel(v: float01) -> float:
    """Colors need to be normalised (from a sRGB space) before computing the relative
    luminance.

    Args:
        v (float): r,g,b channels values between 0 and 1

    Returns:
        float: sRGB channel value for a given rgb color channel
    """
    if v <= 0.04045:
        return v / 12.92
    else:
        return ((v + 0.055) / 1.055) ** 2.4


def relative_luminance(color: RGBColor) -> float:
    """Compute the relative luminance of a color.

    Args:
        color (tuple): rgb color tuple ``(r, g, b)``

    Returns:
        float: relative luminance of a color
    """
    r, g, b = color
    r_ = sRGB_channel(r)
    g_ = sRGB_channel(g)
    b_ = sRGB_channel(b)

    return 0.2126 * r_ + 0.7152 * g_ + 0.0722 * b_


def contrast_ratio(color1: RGBColor, color2: RGBColor) -> float:
    """Compute the contrast ratio between two colors.

    Args:
        color1 (tuple): rgb color tuple ``(r, g, b)``
        color2 (tuple): rgb color tuple ``(r, g, b)``

    Returns:
        float: contrast ratio between two colors
    """

    l1 = relative_luminance(color1)
    l2 = relative_luminance(color2)

    if l1 > l2:
        return (l1 + 0.05) / (l2 + 0.05)
    else:
        return (l2 + 0.05) / (l1 + 0.05)


def passes_contrast(
    color1: RGBColor, color2: RGBColor, level="AA"
) -> Union[bool, float]:
    """Method to verify the contrast ratio between two colours.

    Args:
        color1 (tuple): rgb color tuple ``(r, g, b)``
        color2 (tuple): rgb color tuple ``(r, g, b)``
        level (str, optional): WCAG contrast level. Defaults to "AA".
    """

    ratio = contrast_ratio(color1, color2)

    if level == "AA":
        if ratio >= 4.5:
            return round(ratio, 2)
        else:
            return False
    elif level == "AAA":
        if ratio >= 7.0:
            return round(ratio, 2)
        else:
            return False

    raise ValueError("level must be either 'AA' or 'AAA'")


def get_wcag_level_normal_text(contrast: float) -> str:
    """Does the given contrast meet level AA or level AAA for normal size text"""
    if contrast >= 7:
        return "AAA"
    elif contrast >= 4.5:
        return "AA"
    else:
        return ""


def get_wcag_level_large_text(contrast: float) -> str:
    """Does the given contrast meet level AA or level AAA for large text"""
    if contrast >= 4.5:
        return "AAA"
    elif contrast >= 3:
        return "AA"
    else:
        return ""


def hexstr_without_hash(hex_color: str) -> str:
    """Remove '#' from hex color strings.

    Usage example:

    >>> hexstr_without_hash("#fff")
    'fff'
    """
    return hex_color.replace("#", "")
