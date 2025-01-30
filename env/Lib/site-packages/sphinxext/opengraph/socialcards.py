"""Build a PNG card for each page meant for social media."""
import hashlib
from pathlib import Path
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.image as mpimg
from sphinx.util import logging

matplotlib.use("agg")

LOGGER = logging.getLogger(__name__)
HERE = Path(__file__).parent
MAX_CHAR_PAGE_TITLE = 75
MAX_CHAR_DESCRIPTION = 175

# Default configuration for this functionality
DEFAULT_SOCIAL_CONFIG = {
    "enable": True,
    "site_url": True,
    "site_title": True,
    "page_title": True,
    "description": True,
}


# Default configuration for the figure style
DEFAULT_KWARGS_FIG = {
    "enable": True,
    "site_url": True,
}


# These functions are used when creating social card objects to set MPL values.
# They must be defined here otherwise Sphinx errors when trying to pickle them.
# They are dependent on the `multiple` variable defined when the figure is created.
# Because they are depending on the figure size and renderer used to generate them.
def _set_page_title_line_width():
    return 825


def _set_description_line_width():
    return 1000


def create_social_card(
    app, config_social, site_name, page_title, description, url_text, page_path
):
    """Create a social preview card according to page metadata.

    This uses page metadata and calls a render function to generate the image.
    It also passes configuration through to the rendering function.
    If Matplotlib objects are present in the `app` environment, it reuses them.
    """

    # Add a hash to the image path based on metadata to bust caches
    # ref: https://developer.twitter.com/en/docs/twitter-for-websites/cards/guides/troubleshooting-cards#refreshing_images  # noqa
    hash = hashlib.sha1(
        (site_name + page_title + description + str(config_social)).encode()
    ).hexdigest()[:8]

    # Define the file path we'll use for this image
    path_images_relative = Path("_images/social_previews")
    filename_image = f"summary_{page_path.replace('/', '_')}_{hash}.png"

    # Absolute path used to save the image
    path_images_absolute = Path(app.builder.outdir) / path_images_relative
    path_images_absolute.mkdir(exist_ok=True, parents=True)
    path_image = path_images_absolute / filename_image

    # If the image already exists then we can just skip creating a new one.
    # This is because we hash the values of the text + images in the social card.
    # If the hash doesn't change, it means the output should be the same.
    if path_image.exists():
        return

    # These kwargs are used to generate the base figure image
    kwargs_fig = {}

    # Large image to the top right
    if config_social.get("image"):
        kwargs_fig["image"] = Path(app.builder.srcdir) / config_social.get("image")
    elif app.config.html_logo:
        kwargs_fig["image"] = Path(app.builder.srcdir) / app.config.html_logo

    # Mini image to the bottom right
    if config_social.get("image_mini"):
        kwargs_fig["image_mini"] = Path(app.builder.srcdir) / config_social.get(
            "image_mini"
        )
    else:
        kwargs_fig["image_mini"] = (
            Path(__file__).parent / "_static/sphinx-logo-shadow.png"
        )

    # Validation on the images
    for img in ["image_mini", "image"]:
        impath = kwargs_fig.get(img)
        if not impath:
            continue

        # If image is an SVG replace it with None
        if impath.suffix.lower() == ".svg":
            LOGGER.warning(f"[Social card] %s cannot be an SVG image, skipping...", img)
            kwargs_fig[img] = None

        # If image doesn't exist, throw a warning and replace with none
        if not impath.exists():
            LOGGER.warning(f"[Social card]: %s file doesn't exist, skipping...", img)
            kwargs_fig[img] = None

    # These are passed directly from the user configuration to our plotting function
    pass_through_config = ["text_color", "line_color", "background_color", "font"]
    for config in pass_through_config:
        if config_social.get(config):
            kwargs_fig[config] = config_social.get(config)

    # Generate the image and store the matplotlib objects so that we can re-use them
    if hasattr(app.env, "ogp_social_card_plt_objects"):
        plt_objects = app.env.ogp_social_card_plt_objects
    else:
        plt_objects = None
    plt_objects = render_social_card(
        path_image,
        site_name,
        page_title,
        description,
        url_text,
        plt_objects,
        kwargs_fig,
    )
    app.env.ogp_social_card_plt_objects = plt_objects

    # Path relative to build folder will be what we use for linking the URL
    path_relative_to_build = path_images_relative / filename_image
    return path_relative_to_build


def render_social_card(
    path,
    site_title=None,
    page_title=None,
    description=None,
    siteurl=None,
    plt_objects=None,
    kwargs_fig=None,
):
    """Render a social preview card with Matplotlib and write to disk."""
    # If objects is None it means this is the first time plotting.
    # Create the figure objects and return them so that we re-use them later.
    if plt_objects is None:
        (
            fig,
            txt_site_title,
            txt_page_title,
            txt_description,
            txt_url,
        ) = create_social_card_objects(**kwargs_fig)
    else:
        fig, txt_site_title, txt_page_title, txt_description, txt_url = plt_objects

    # Update the matplotlib text objects with new text from this page
    txt_site_title.set_text(site_title)
    txt_page_title.set_text(page_title)
    txt_description.set_text(description)
    txt_url.set_text(siteurl)

    # Save the image
    fig.savefig(path, facecolor=None)
    return fig, txt_site_title, txt_page_title, txt_description, txt_url


def create_social_card_objects(
    image=None,
    image_mini=None,
    page_title_color="#2f363d",
    description_color="#585e63",
    site_title_color="#585e63",
    site_url_color="#2f363d",
    background_color="white",
    line_color="#5A626B",
    font=None,
):
    """Create the Matplotlib objects for the first time."""
    # If no font specified, load the Roboto Flex font as a fallback
    if font is None:
        path_font = Path(__file__).parent / "_static/Roboto-Flex.ttf"
        roboto_font = matplotlib.font_manager.FontEntry(
            fname=str(path_font), name="Roboto Flex"
        )
        matplotlib.font_manager.fontManager.addfont(path_font)
        font = roboto_font.name

    # Because Matplotlib doesn't let you specify figures in pixels, only inches
    # This `multiple` results in a scale of about 1146px by 600px
    # Which is roughly the recommended size for OpenGraph images
    # ref: https://opengraph.xyz
    ratio = 1200 / 628
    multiple = 6
    fig = plt.figure(figsize=(ratio * multiple, multiple))
    fig.set_facecolor(background_color)

    # Text axis
    axtext = fig.add_axes((0, 0, 1, 1))

    # Image axis
    ax_x, ax_y, ax_w, ax_h = (0.65, 0.65, 0.3, 0.3)
    axim_logo = fig.add_axes((ax_x, ax_y, ax_w, ax_h), anchor="NE")

    # Image mini axis
    ax_x, ax_y, ax_w, ax_h = (0.82, 0.1, 0.1, 0.1)
    axim_mini = fig.add_axes((ax_x, ax_y, ax_w, ax_h), anchor="NE")

    # Line at the bottom axis
    axline = fig.add_axes((-0.1, -0.04, 1.2, 0.1))

    # Axes configuration
    left_margin = 0.05
    with plt.rc_context({"font.family": font}):
        # Site title
        # Smaller font, just above page title
        site_title_y_offset = 0.87
        txt_site = axtext.text(
            left_margin,
            site_title_y_offset,
            "Test site title",
            {
                "size": 24,
            },
            ha="left",
            va="top",
            wrap=True,
            c=site_title_color,
        )

        # Page title
        # A larger font for more visibility
        page_title_y_offset = 0.77

        txt_page = axtext.text(
            left_margin,
            page_title_y_offset,
            "Test page title, a bit longer to demo",
            {"size": 46, "color": "k", "fontweight": "bold"},
            ha="left",
            va="top",
            wrap=True,
            c=page_title_color,
        )

        txt_page._get_wrap_line_width = _set_page_title_line_width

        # description
        # Just below site title, smallest font and many lines.
        # Our target length is 160 characters, so it should be
        # two lines at full width with some room to spare at this length.
        description_y_offset = 0.2
        txt_description = axtext.text(
            left_margin,
            description_y_offset,
            (
                "A longer description that we use to ,"
                "show off what the descriptions look like."
            ),
            {"size": 17},
            ha="left",
            va="bottom",
            wrap=True,
            c=description_color,
        )
        txt_description._get_wrap_line_width = _set_description_line_width

        # url
        # Aligned to the left of the mini image
        url_y_axis_ofset = 0.12
        txt_url = axtext.text(
            left_margin,
            url_y_axis_ofset,
            "testurl.org",
            {"size": 22},
            ha="left",
            va="bottom",
            fontweight="bold",
            c=site_url_color,
        )

    if isinstance(image_mini, Path):
        img = mpimg.imread(image_mini)
        axim_mini.imshow(img)

    # Put the logo in the top right if it exists
    if isinstance(image, Path):
        img = mpimg.imread(image)
        yw, xw = img.shape[:2]

        # Axis is square and width is longest image axis
        longest = max([yw, xw])
        axim_logo.set_xlim([0, longest])
        axim_logo.set_ylim([longest, 0])

        # Center it on the non-long axis
        xdiff = (longest - xw) / 2
        ydiff = (longest - yw) / 2
        axim_logo.imshow(img, extent=[xdiff, xw + xdiff, yw + ydiff, ydiff])

    # Put a colored line at the bottom of the figure
    axline.hlines(0, 0, 1, lw=25, color=line_color)

    # Remove the ticks and borders from all axes for a clean look
    for ax in fig.axes:
        ax.set_axis_off()
    return fig, txt_site, txt_page, txt_description, txt_url
