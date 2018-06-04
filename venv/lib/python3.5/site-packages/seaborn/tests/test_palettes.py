import warnings
import colorsys
import numpy as np
import matplotlib as mpl

import nose.tools as nt
import numpy.testing as npt

from .. import palettes, utils, rcmod
from ..external import husl
from ..xkcd_rgb import xkcd_rgb
from ..crayons import crayons


class TestColorPalettes(object):

    def test_current_palette(self):

        pal = palettes.color_palette(["red", "blue", "green"], 3)
        rcmod.set_palette(pal, 3)
        nt.assert_equal(pal, utils.get_color_cycle())
        rcmod.set()

    def test_palette_context(self):

        default_pal = palettes.color_palette()
        context_pal = palettes.color_palette("muted")

        with palettes.color_palette(context_pal):
            nt.assert_equal(utils.get_color_cycle(), context_pal)

        nt.assert_equal(utils.get_color_cycle(), default_pal)

    def test_big_palette_context(self):

        original_pal = palettes.color_palette("deep", n_colors=8)
        context_pal = palettes.color_palette("husl", 10)

        rcmod.set_palette(original_pal)
        with palettes.color_palette(context_pal, 10):
            nt.assert_equal(utils.get_color_cycle(), context_pal)

        nt.assert_equal(utils.get_color_cycle(), original_pal)

        # Reset default
        rcmod.set()

    def test_seaborn_palettes(self):

        pals = "deep", "muted", "pastel", "bright", "dark", "colorblind"
        for name in pals:
            pal_out = palettes.color_palette(name)
            nt.assert_equal(len(pal_out), 6)

    def test_hls_palette(self):

        hls_pal1 = palettes.hls_palette()
        hls_pal2 = palettes.color_palette("hls")
        npt.assert_array_equal(hls_pal1, hls_pal2)

    def test_husl_palette(self):

        husl_pal1 = palettes.husl_palette()
        husl_pal2 = palettes.color_palette("husl")
        npt.assert_array_equal(husl_pal1, husl_pal2)

    def test_mpl_palette(self):

        mpl_pal1 = palettes.mpl_palette("Reds")
        mpl_pal2 = palettes.color_palette("Reds")
        npt.assert_array_equal(mpl_pal1, mpl_pal2)

    def test_mpl_dark_palette(self):

        mpl_pal1 = palettes.mpl_palette("Blues_d")
        mpl_pal2 = palettes.color_palette("Blues_d")
        npt.assert_array_equal(mpl_pal1, mpl_pal2)

    def test_bad_palette_name(self):

        with nt.assert_raises(ValueError):
            palettes.color_palette("IAmNotAPalette")

    def test_terrible_palette_name(self):

        with nt.assert_raises(ValueError):
            palettes.color_palette("jet")

    def test_bad_palette_colors(self):

        pal = ["red", "blue", "iamnotacolor"]
        with nt.assert_raises(ValueError):
            palettes.color_palette(pal)

    def test_palette_desat(self):

        pal1 = palettes.husl_palette(6)
        pal1 = [utils.desaturate(c, .5) for c in pal1]
        pal2 = palettes.color_palette("husl", desat=.5)
        npt.assert_array_equal(pal1, pal2)

    def test_palette_is_list_of_tuples(self):

        pal_in = np.array(["red", "blue", "green"])
        pal_out = palettes.color_palette(pal_in, 3)

        nt.assert_is_instance(pal_out, list)
        nt.assert_is_instance(pal_out[0], tuple)
        nt.assert_is_instance(pal_out[0][0], float)
        nt.assert_equal(len(pal_out[0]), 3)

    def test_palette_cycles(self):

        deep = palettes.color_palette("deep")
        double_deep = palettes.color_palette("deep", 12)
        nt.assert_equal(double_deep, deep + deep)

    def test_hls_values(self):

        pal1 = palettes.hls_palette(6, h=0)
        pal2 = palettes.hls_palette(6, h=.5)
        pal2 = pal2[3:] + pal2[:3]
        npt.assert_array_almost_equal(pal1, pal2)

        pal_dark = palettes.hls_palette(5, l=.2)
        pal_bright = palettes.hls_palette(5, l=.8)
        npt.assert_array_less(list(map(sum, pal_dark)),
                              list(map(sum, pal_bright)))

        pal_flat = palettes.hls_palette(5, s=.1)
        pal_bold = palettes.hls_palette(5, s=.9)
        npt.assert_array_less(list(map(np.std, pal_flat)),
                              list(map(np.std, pal_bold)))

    def test_husl_values(self):

        pal1 = palettes.husl_palette(6, h=0)
        pal2 = palettes.husl_palette(6, h=.5)
        pal2 = pal2[3:] + pal2[:3]
        npt.assert_array_almost_equal(pal1, pal2)

        pal_dark = palettes.husl_palette(5, l=.2)
        pal_bright = palettes.husl_palette(5, l=.8)
        npt.assert_array_less(list(map(sum, pal_dark)),
                              list(map(sum, pal_bright)))

        pal_flat = palettes.husl_palette(5, s=.1)
        pal_bold = palettes.husl_palette(5, s=.9)
        npt.assert_array_less(list(map(np.std, pal_flat)),
                              list(map(np.std, pal_bold)))

    def test_cbrewer_qual(self):

        pal_short = palettes.mpl_palette("Set1", 4)
        pal_long = palettes.mpl_palette("Set1", 6)
        nt.assert_equal(pal_short, pal_long[:4])

        pal_full = palettes.mpl_palette("Set2", 8)
        pal_long = palettes.mpl_palette("Set2", 10)
        nt.assert_equal(pal_full, pal_long[:8])

    def test_mpl_reversal(self):

        pal_forward = palettes.mpl_palette("BuPu", 6)
        pal_reverse = palettes.mpl_palette("BuPu_r", 6)
        npt.assert_array_almost_equal(pal_forward, pal_reverse[::-1])

    def test_rgb_from_hls(self):

        color = .5, .8, .4
        rgb_got = palettes._color_to_rgb(color, "hls")
        rgb_want = colorsys.hls_to_rgb(*color)
        nt.assert_equal(rgb_got, rgb_want)

    def test_rgb_from_husl(self):

        color = 120, 50, 40
        rgb_got = palettes._color_to_rgb(color, "husl")
        rgb_want = husl.husl_to_rgb(*color)
        nt.assert_equal(rgb_got, rgb_want)

    def test_rgb_from_xkcd(self):

        color = "dull red"
        rgb_got = palettes._color_to_rgb(color, "xkcd")
        rgb_want = xkcd_rgb[color]
        nt.assert_equal(rgb_got, rgb_want)

    def test_light_palette(self):

        pal_forward = palettes.light_palette("red")
        pal_reverse = palettes.light_palette("red", reverse=True)
        npt.assert_array_almost_equal(pal_forward, pal_reverse[::-1])

        red = tuple(mpl.colors.colorConverter.to_rgba("red"))
        nt.assert_equal(tuple(pal_forward[-1]), red)

        pal_cmap = palettes.light_palette("blue", as_cmap=True)
        nt.assert_is_instance(pal_cmap, mpl.colors.LinearSegmentedColormap)

    def test_dark_palette(self):

        pal_forward = palettes.dark_palette("red")
        pal_reverse = palettes.dark_palette("red", reverse=True)
        npt.assert_array_almost_equal(pal_forward, pal_reverse[::-1])

        red = tuple(mpl.colors.colorConverter.to_rgba("red"))
        nt.assert_equal(tuple(pal_forward[-1]), red)

        pal_cmap = palettes.dark_palette("blue", as_cmap=True)
        nt.assert_is_instance(pal_cmap, mpl.colors.LinearSegmentedColormap)

    def test_blend_palette(self):

        colors = ["red", "yellow", "white"]
        pal_cmap = palettes.blend_palette(colors, as_cmap=True)
        nt.assert_is_instance(pal_cmap, mpl.colors.LinearSegmentedColormap)

    def test_cubehelix_against_matplotlib(self):

        x = np.linspace(0, 1, 8)
        mpl_pal = mpl.cm.cubehelix(x)[:, :3].tolist()

        sns_pal = palettes.cubehelix_palette(8, start=0.5, rot=-1.5, hue=1,
                                             dark=0, light=1, reverse=True)

        nt.assert_list_equal(sns_pal, mpl_pal)

    def test_cubehelix_n_colors(self):

        for n in [3, 5, 8]:
            pal = palettes.cubehelix_palette(n)
            nt.assert_equal(len(pal), n)

    def test_cubehelix_reverse(self):

        pal_forward = palettes.cubehelix_palette()
        pal_reverse = palettes.cubehelix_palette(reverse=True)
        nt.assert_list_equal(pal_forward, pal_reverse[::-1])

    def test_cubehelix_cmap(self):

        cmap = palettes.cubehelix_palette(as_cmap=True)
        nt.assert_is_instance(cmap, mpl.colors.ListedColormap)
        pal = palettes.cubehelix_palette()
        x = np.linspace(0, 1, 6)
        npt.assert_array_equal(cmap(x)[:, :3], pal)

        cmap_rev = palettes.cubehelix_palette(as_cmap=True, reverse=True)
        x = np.linspace(0, 1, 6)
        pal_forward = cmap(x).tolist()
        pal_reverse = cmap_rev(x[::-1]).tolist()
        nt.assert_list_equal(pal_forward, pal_reverse)

    def test_xkcd_palette(self):

        names = list(xkcd_rgb.keys())[10:15]
        colors = palettes.xkcd_palette(names)
        for name, color in zip(names, colors):
            as_hex = mpl.colors.rgb2hex(color)
            nt.assert_equal(as_hex, xkcd_rgb[name])

    def test_crayon_palette(self):

        names = list(crayons.keys())[10:15]
        colors = palettes.crayon_palette(names)
        for name, color in zip(names, colors):
            as_hex = mpl.colors.rgb2hex(color)
            nt.assert_equal(as_hex, crayons[name].lower())

    def test_color_codes(self):

        palettes.set_color_codes("deep")
        colors = palettes.color_palette("deep") + [".1"]
        for code, color in zip("bgrmyck", colors):
            rgb_want = mpl.colors.colorConverter.to_rgb(color)
            rgb_got = mpl.colors.colorConverter.to_rgb(code)
            nt.assert_equal(rgb_want, rgb_got)
        palettes.set_color_codes("reset")

    def test_as_hex(self):

        pal = palettes.color_palette("deep")
        for rgb, hex in zip(pal, pal.as_hex()):
            nt.assert_equal(mpl.colors.rgb2hex(rgb), hex)

    def test_preserved_palette_length(self):

        pal_in = palettes.color_palette("Set1", 10)
        pal_out = palettes.color_palette(pal_in)
        nt.assert_equal(pal_in, pal_out)

    def test_get_color_cycle(self):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            result = utils.get_color_cycle()
            expected = mpl.rcParams['axes.color_cycle']
        nt.assert_equal(result, expected)
