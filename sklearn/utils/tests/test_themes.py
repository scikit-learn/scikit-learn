import pytest

from sklearn.utils import themes

TEMPLATE_SIMPLE = """
#$id {
  color: $$color_1;
  background-color: $$color_2;
}
""".replace("  ", "").replace("\n", "")  # noqa

TEMPLATE_SIMPLE_2 = """
#$id {
  color: $$color_1;
}
#$id div {
  background-color: $$color_2;
}
""".replace("  ", "").replace("\n", "")  # noqa


@pytest.mark.parametrize(
    "color_1,color_2,expected",
    [
        ("red", "#131313", "#$id {color: red;background-color: #131313;}"),
        ("blue", "#070707", "#$id {color: blue;background-color: #070707;}"),
    ],
)
def test_css_template_simple(color_1, color_2, expected):
    # Test CSS Template rendering.

    tmpl = themes.CssTemplate(TEMPLATE_SIMPLE)
    result = tmpl.substitute(color_1=color_1, color_2=color_2)

    assert expected == result


@pytest.mark.parametrize(
    "tmpl,expected",
    [
        (TEMPLATE_SIMPLE, "#$id {color: red;background-color: blue;}"),
        (TEMPLATE_SIMPLE_2, "#$id {color: red;}#$id div {background-color: blue;}"),
    ],
)
def test_theme_builder(tmpl, expected):
    # Test theme builder renders correctly.

    result = themes.theme_builder(
        template=themes.CssTemplate(tmpl), color_1="red", color_2="blue"
    )

    assert expected == result
