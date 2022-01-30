import colorsys
import numpy as np
import matplotlib as mpl

import pytest
import numpy.testing as npt

from .. import palettes, utils, rcmod
from ..external import husl
from ..colors import xkcd_rgb, crayons


class TestColorPalettes:

    def test_current_palette(self):

        pal = palettes.color_palette(["red", "blue", "green"])
        rcmod.set_palette(pal)
        assert pal == utils.get_color_cycle()
        rcmod.set()

    def test_palette_context(self):

        default_pal = palettes.color_palette()
        context_pal = palettes.color_palette("muted")

        with palettes.color_palette(context_pal):
            assert utils.get_color_cycle() == context_pal

        assert utils.get_color_cycle() == default_pal

    def test_big_palette_context(self):

        original_pal = palettes.color_palette("deep", n_colors=8)
        context_pal = palettes.color_palette("husl", 10)

        rcmod.set_palette(original_pal)
        with palettes.color_palette(context_pal, 10):
            assert utils.get_color_cycle() == context_pal

        assert utils.get_color_cycle() == original_pal

        # Reset default
        rcmod.set()

    def test_palette_size(self):

        pal = palettes.color_palette("deep")
        assert len(pal) == palettes.QUAL_PALETTE_SIZES["deep"]

        pal = palettes.color_palette("pastel6")
        assert len(pal) == palettes.QUAL_PALETTE_SIZES["pastel6"]

        pal = palettes.color_palette("Set3")
        assert len(pal) == palettes.QUAL_PALETTE_SIZES["Set3"]

        pal = palettes.color_palette("husl")
        assert len(pal) == 6

        pal = palettes.color_palette("Greens")
        assert len(pal) == 6

    def test_seaborn_palettes(self):

        pals = "deep", "muted", "pastel", "bright", "dark", "colorblind"
        for name in pals:
            full = palettes.color_palette(name, 10).as_hex()
            short = palettes.color_palette(name + "6", 6).as_hex()
            b, _, g, r, m, _, _, _, y, c = full
            assert [b, g, r, m, y, c] == list(short)

    def test_hls_palette(self):

        pal1 = palettes.hls_palette()
        pal2 = palettes.color_palette("hls")
        npt.assert_array_equal(pal1, pal2)

        cmap1 = palettes.hls_palette(as_cmap=True)
        cmap2 = palettes.color_palette("hls", as_cmap=True)
        npt.assert_array_equal(cmap1([.2, .8]), cmap2([.2, .8]))

    def test_husl_palette(self):

        pal1 = palettes.husl_palette()
        pal2 = palettes.color_palette("husl")
        npt.assert_array_equal(pal1, pal2)

        cmap1 = palettes.husl_palette(as_cmap=True)
        cmap2 = palettes.color_palette("husl", as_cmap=True)
        npt.assert_array_equal(cmap1([.2, .8]), cmap2([.2, .8]))

    def test_mpl_palette(self):

        pal1 = palettes.mpl_palette("Reds")
        pal2 = palettes.color_palette("Reds")
        npt.assert_array_equal(pal1, pal2)

        cmap1 = mpl.cm.get_cmap("Reds")
        cmap2 = palettes.mpl_palette("Reds", as_cmap=True)
        cmap3 = palettes.color_palette("Reds", as_cmap=True)
        npt.assert_array_equal(cmap1, cmap2)
        npt.assert_array_equal(cmap1, cmap3)

    def test_mpl_dark_palette(self):

        mpl_pal1 = palettes.mpl_palette("Blues_d")
        mpl_pal2 = palettes.color_palette("Blues_d")
        npt.assert_array_equal(mpl_pal1, mpl_pal2)

        mpl_pal1 = palettes.mpl_palette("Blues_r_d")
        mpl_pal2 = palettes.color_palette("Blues_r_d")
        npt.assert_array_equal(mpl_pal1, mpl_pal2)

    def test_bad_palette_name(self):

        with pytest.raises(ValueError):
            palettes.color_palette("IAmNotAPalette")

    def test_terrible_palette_name(self):

        with pytest.raises(ValueError):
            palettes.color_palette("jet")

    def test_bad_palette_colors(self):

        pal = ["red", "blue", "iamnotacolor"]
        with pytest.raises(ValueError):
            palettes.color_palette(pal)

    def test_palette_desat(self):

        pal1 = palettes.husl_palette(6)
        pal1 = [utils.desaturate(c, .5) for c in pal1]
        pal2 = palettes.color_palette("husl", desat=.5)
        npt.assert_array_equal(pal1, pal2)

    def test_palette_is_list_of_tuples(self):

        pal_in = np.array(["red", "blue", "green"])
        pal_out = palettes.color_palette(pal_in, 3)

        assert isinstance(pal_out, list)
        assert isinstance(pal_out[0], tuple)
        assert isinstance(pal_out[0][0], float)
        assert len(pal_out[0]) == 3

    def test_palette_cycles(self):

        deep = palettes.color_palette("deep6")
        double_deep = palettes.color_palette("deep6", 12)
        assert double_deep == deep + deep

    def test_hls_values(self):

        pal1 = palettes.hls_palette(6, h=0)
        pal2 = palettes.hls_palette(6, h=.5)
        pal2 = pal2[3:] + pal2[:3]
        npt.assert_array_almost_equal(pal1, pal2)

        pal_dark = palettes.hls_palette(5, l=.2)  # noqa
        pal_bright = palettes.hls_palette(5, l=.8)  # noqa
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

        pal_dark = palettes.husl_palette(5, l=.2)  # noqa
        pal_bright = palettes.husl_palette(5, l=.8)  # noqa
        npt.assert_array_less(list(map(sum, pal_dark)),
                              list(map(sum, pal_bright)))

        pal_flat = palettes.husl_palette(5, s=.1)
        pal_bold = palettes.husl_palette(5, s=.9)
        npt.assert_array_less(list(map(np.std, pal_flat)),
                              list(map(np.std, pal_bold)))

    def test_cbrewer_qual(self):

        pal_short = palettes.mpl_palette("Set1", 4)
        pal_long = palettes.mpl_palette("Set1", 6)
        assert pal_short == pal_long[:4]

        pal_full = palettes.mpl_palette("Set2", 8)
        pal_long = palettes.mpl_palette("Set2", 10)
        assert pal_full == pal_long[:8]

    def test_mpl_reversal(self):

        pal_forward = palettes.mpl_palette("BuPu", 6)
        pal_reverse = palettes.mpl_palette("BuPu_r", 6)
        npt.assert_array_almost_equal(pal_forward, pal_reverse[::-1])

    def test_rgb_from_hls(self):

        color = .5, .8, .4
        rgb_got = palettes._color_to_rgb(color, "hls")
        rgb_want = colorsys.hls_to_rgb(*color)
        assert rgb_got == rgb_want

    def test_rgb_from_husl(self):

        color = 120, 50, 40
        rgb_got = palettes._color_to_rgb(color, "husl")
        rgb_want = tuple(husl.husl_to_rgb(*color))
        assert rgb_got == rgb_want

        for h in range(0, 360):
            color = h, 100, 100
            rgb = palettes._color_to_rgb(color, "husl")
            assert min(rgb) >= 0
            assert max(rgb) <= 1

    def test_rgb_from_xkcd(self):

        color = "dull red"
        rgb_got = palettes._color_to_rgb(color, "xkcd")
        rgb_want = mpl.colors.to_rgb(xkcd_rgb[color])
        assert rgb_got == rgb_want

    def test_light_palette(self):

        n = 4
        pal_forward = palettes.light_palette("red", n)
        pal_reverse = palettes.light_palette("red", n, reverse=True)
        assert np.allclose(pal_forward, pal_reverse[::-1])

        red = mpl.colors.colorConverter.to_rgb("red")
        assert pal_forward[-1] == red

        pal_f_from_string = palettes.color_palette("light:red", n)
        assert pal_forward[3] == pal_f_from_string[3]

        pal_r_from_string = palettes.color_palette("light:red_r", n)
        assert pal_reverse[3] == pal_r_from_string[3]

        pal_cmap = palettes.light_palette("blue", as_cmap=True)
        assert isinstance(pal_cmap, mpl.colors.LinearSegmentedColormap)

        pal_cmap_from_string = palettes.color_palette("light:blue", as_cmap=True)
        assert pal_cmap(.8) == pal_cmap_from_string(.8)

        pal_cmap = palettes.light_palette("blue", as_cmap=True, reverse=True)
        pal_cmap_from_string = palettes.color_palette("light:blue_r", as_cmap=True)
        assert pal_cmap(.8) == pal_cmap_from_string(.8)

    def test_dark_palette(self):

        n = 4
        pal_forward = palettes.dark_palette("red", n)
        pal_reverse = palettes.dark_palette("red", n, reverse=True)
        assert np.allclose(pal_forward, pal_reverse[::-1])

        red = mpl.colors.colorConverter.to_rgb("red")
        assert pal_forward[-1] == red

        pal_f_from_string = palettes.color_palette("dark:red", n)
        assert pal_forward[3] == pal_f_from_string[3]

        pal_r_from_string = palettes.color_palette("dark:red_r", n)
        assert pal_reverse[3] == pal_r_from_string[3]

        pal_cmap = palettes.dark_palette("blue", as_cmap=True)
        assert isinstance(pal_cmap, mpl.colors.LinearSegmentedColormap)

        pal_cmap_from_string = palettes.color_palette("dark:blue", as_cmap=True)
        assert pal_cmap(.8) == pal_cmap_from_string(.8)

        pal_cmap = palettes.dark_palette("blue", as_cmap=True, reverse=True)
        pal_cmap_from_string = palettes.color_palette("dark:blue_r", as_cmap=True)
        assert pal_cmap(.8) == pal_cmap_from_string(.8)

    def test_diverging_palette(self):

        h_neg, h_pos = 100, 200
        sat, lum = 70, 50
        args = h_neg, h_pos, sat, lum

        n = 12
        pal = palettes.diverging_palette(*args, n=n)
        neg_pal = palettes.light_palette((h_neg, sat, lum), int(n // 2),
                                         input="husl")
        pos_pal = palettes.light_palette((h_pos, sat, lum), int(n // 2),
                                         input="husl")
        assert len(pal) == n
        assert pal[0] == neg_pal[-1]
        assert pal[-1] == pos_pal[-1]

        pal_dark = palettes.diverging_palette(*args, n=n, center="dark")
        assert np.mean(pal[int(n / 2)]) > np.mean(pal_dark[int(n / 2)])

        pal_cmap = palettes.diverging_palette(*args, as_cmap=True)
        assert isinstance(pal_cmap, mpl.colors.LinearSegmentedColormap)

    def test_blend_palette(self):

        colors = ["red", "yellow", "white"]
        pal_cmap = palettes.blend_palette(colors, as_cmap=True)
        assert isinstance(pal_cmap, mpl.colors.LinearSegmentedColormap)

        colors = ["red", "blue"]
        pal = palettes.blend_palette(colors)
        pal_str = "blend:" + ",".join(colors)
        pal_from_str = palettes.color_palette(pal_str)
        assert pal == pal_from_str

    def test_cubehelix_against_matplotlib(self):

        x = np.linspace(0, 1, 8)
        mpl_pal = mpl.cm.cubehelix(x)[:, :3].tolist()

        sns_pal = palettes.cubehelix_palette(8, start=0.5, rot=-1.5, hue=1,
                                             dark=0, light=1, reverse=True)

        assert sns_pal == mpl_pal

    def test_cubehelix_n_colors(self):

        for n in [3, 5, 8]:
            pal = palettes.cubehelix_palette(n)
            assert len(pal) == n

    def test_cubehelix_reverse(self):

        pal_forward = palettes.cubehelix_palette()
        pal_reverse = palettes.cubehelix_palette(reverse=True)
        assert pal_forward == pal_reverse[::-1]

    def test_cubehelix_cmap(self):

        cmap = palettes.cubehelix_palette(as_cmap=True)
        assert isinstance(cmap, mpl.colors.ListedColormap)
        pal = palettes.cubehelix_palette()
        x = np.linspace(0, 1, 6)
        npt.assert_array_equal(cmap(x)[:, :3], pal)

        cmap_rev = palettes.cubehelix_palette(as_cmap=True, reverse=True)
        x = np.linspace(0, 1, 6)
        pal_forward = cmap(x).tolist()
        pal_reverse = cmap_rev(x[::-1]).tolist()
        assert pal_forward == pal_reverse

    def test_cubehelix_code(self):

        color_palette = palettes.color_palette
        cubehelix_palette = palettes.cubehelix_palette

        pal1 = color_palette("ch:", 8)
        pal2 = color_palette(cubehelix_palette(8))
        assert pal1 == pal2

        pal1 = color_palette("ch:.5, -.25,hue = .5,light=.75", 8)
        pal2 = color_palette(cubehelix_palette(8, .5, -.25, hue=.5, light=.75))
        assert pal1 == pal2

        pal1 = color_palette("ch:h=1,r=.5", 9)
        pal2 = color_palette(cubehelix_palette(9, hue=1, rot=.5))
        assert pal1 == pal2

        pal1 = color_palette("ch:_r", 6)
        pal2 = color_palette(cubehelix_palette(6, reverse=True))
        assert pal1 == pal2

        pal1 = color_palette("ch:_r", as_cmap=True)
        pal2 = cubehelix_palette(6, reverse=True, as_cmap=True)
        assert pal1(.5) == pal2(.5)

    def test_xkcd_palette(self):

        names = list(xkcd_rgb.keys())[10:15]
        colors = palettes.xkcd_palette(names)
        for name, color in zip(names, colors):
            as_hex = mpl.colors.rgb2hex(color)
            assert as_hex == xkcd_rgb[name]

    def test_crayon_palette(self):

        names = list(crayons.keys())[10:15]
        colors = palettes.crayon_palette(names)
        for name, color in zip(names, colors):
            as_hex = mpl.colors.rgb2hex(color)
            assert as_hex == crayons[name].lower()

    def test_color_codes(self):

        palettes.set_color_codes("deep")
        colors = palettes.color_palette("deep6") + [".1"]
        for code, color in zip("bgrmyck", colors):
            rgb_want = mpl.colors.colorConverter.to_rgb(color)
            rgb_got = mpl.colors.colorConverter.to_rgb(code)
            assert rgb_want == rgb_got
        palettes.set_color_codes("reset")

        with pytest.raises(ValueError):
            palettes.set_color_codes("Set1")

    def test_as_hex(self):

        pal = palettes.color_palette("deep")
        for rgb, hex in zip(pal, pal.as_hex()):
            assert mpl.colors.rgb2hex(rgb) == hex

    def test_preserved_palette_length(self):

        pal_in = palettes.color_palette("Set1", 10)
        pal_out = palettes.color_palette(pal_in)
        assert pal_in == pal_out

    def test_html_rep(self):

        pal = palettes.color_palette()
        html = pal._repr_html_()
        for color in pal.as_hex():
            assert color in html
