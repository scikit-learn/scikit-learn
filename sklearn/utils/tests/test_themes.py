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


def test_original_style():
    # Test that themes.LIGHT matches up with the original stylesheet.
    # This test can be deleted when themes.LIGHT style changes.
    _STYLE = """
    #$id {
      color: black;
    }
    #$id pre{
      padding: 0;
    }
    #$id div.sk-toggleable {
      background-color: white;
    }
    #$id label.sk-toggleable__label {
      cursor: pointer;
      display: block;
      width: 100%;
      margin-bottom: 0;
      padding: 0.3em;
      box-sizing: border-box;
      text-align: center;
    }
    #$id label.sk-toggleable__label-arrow:before {
      content: "▸";
      float: left;
      margin-right: 0.25em;
      color: #696969;
    }
    #$id label.sk-toggleable__label-arrow:hover:before {
      color: black;
    }
    #$id div.sk-estimator:hover label.sk-toggleable__label-arrow:before {
      color: black;
    }
    #$id div.sk-toggleable__content {
      max-height: 0;
      max-width: 0;
      overflow: hidden;
      text-align: left;
      background-color: #f0f8ff;
    }
    #$id div.sk-toggleable__content pre {
      margin: 0.2em;
      color: black;
      border-radius: 0.25em;
      background-color: #f0f8ff;
    }
    #$id input.sk-toggleable__control:checked~div.sk-toggleable__content {
      max-height: 200px;
      max-width: 100%;
      overflow: auto;
    }
    #$id input.sk-toggleable__control:checked~label.sk-toggleable__label-arrow:before {
      content: "▾";
    }
    #$id div.sk-estimator input.sk-toggleable__control:checked~label.sk-toggleable__label {
      background-color: #d4ebff;
    }
    #$id div.sk-label input.sk-toggleable__control:checked~label.sk-toggleable__label {
      background-color: #d4ebff;
    }
    #$id input.sk-hidden--visually {
      border: 0;
      clip: rect(1px 1px 1px 1px);
      clip: rect(1px, 1px, 1px, 1px);
      height: 1px;
      margin: -1px;
      overflow: hidden;
      padding: 0;
      position: absolute;
      width: 1px;
    }
    #$id div.sk-estimator {
      font-family: monospace;
      background-color: #f0f8ff;
      border: 1px dotted black;
      border-radius: 0.25em;
      box-sizing: border-box;
      margin-bottom: 0.5em;
    }
    #$id div.sk-estimator:hover {
      background-color: #d4ebff;
    }
    #$id div.sk-parallel-item::after {
      content: "";
      width: 100%;
      border-bottom: 1px solid gray;
      flex-grow: 1;
    }
    #$id div.sk-label:hover label.sk-toggleable__label {
      background-color: #d4ebff;
    }
    #$id div.sk-serial::before {
      content: "";
      position: absolute;
      border-left: 1px solid gray;
      box-sizing: border-box;
      top: 0;
      bottom: 0;
      left: 50%;
      z-index: 0;
    }
    #$id div.sk-serial {
      display: flex;
      flex-direction: column;
      align-items: center;
      background-color: white;
      padding-right: 0.2em;
      padding-left: 0.2em;
      position: relative;
    }
    #$id div.sk-item {
      position: relative;
      z-index: 1;
    }
    #$id div.sk-parallel {
      display: flex;
      align-items: stretch;
      justify-content: center;
      background-color: white;
      position: relative;
    }
    #$id div.sk-item::before, #$id div.sk-parallel-item::before {
      content: "";
      position: absolute;
      border-left: 1px solid gray;
      box-sizing: border-box;
      top: 0;
      bottom: 0;
      left: 50%;
      z-index: -1;
    }
    #$id div.sk-parallel-item {
      display: flex;
      flex-direction: column;
      z-index: 1;
      position: relative;
      background-color: white;
    }
    #$id div.sk-parallel-item:first-child::after {
      align-self: flex-end;
      width: 50%;
    }
    #$id div.sk-parallel-item:last-child::after {
      align-self: flex-start;
      width: 50%;
    }
    #$id div.sk-parallel-item:only-child::after {
      width: 0;
    }
    #$id div.sk-dashed-wrapped {
      border: 1px dashed gray;
      margin: 0 0.4em 0.5em 0.4em;
      box-sizing: border-box;
      padding-bottom: 0.4em;
      background-color: white;
    }
    #$id div.sk-label label {
      font-family: monospace;
      font-weight: bold;
      display: inline-block;
      line-height: 1.2em;
    }
    #$id div.sk-label-container {
      text-align: center;
    }
    #$id div.sk-container {
      /* jupyter's `normalize.less` sets `[hidden] { display: none; }`
         but bootstrap.min.css set `[hidden] { display: none !important; }`
         so we also need the `!important` here to be able to override the
         default hidden behavior on the sphinx rendered scikit-learn.org.
         See: https://github.com/scikit-learn/scikit-learn/issues/21755 */
      display: inline-block !important;
      position: relative;
    }
    #$id div.sk-text-repr-fallback {
      display: none;
    }
    """.replace("  ", "").replace("\n", "")  # noqa

    assert themes.LIGHT == _STYLE


@pytest.mark.parametrize("color_1,color_2", [("red", "#131313"), ("blue", "#070707")])
def test_css_template_simple(color_1, color_2):
    # Test CSS Template rendering.

    tmpl = themes.CssTemplate(TEMPLATE_SIMPLE)
    result = tmpl.substitute(color_1=color_1, color_2=color_2)

    if color_1 == "red":
        assert "#$id {color: red;background-color: #131313;}" == result
    elif color_1 == "blue":
        assert "#$id {color: blue;background-color: #070707;}" == result


@pytest.mark.parametrize("tmpl", [TEMPLATE_SIMPLE, TEMPLATE_SIMPLE_2])
def test_theme_builder(tmpl):
    # Test theme builder renders correctly.

    result = themes.theme_builder(
        template=themes.CssTemplate(tmpl), color_1="red", color_2="blue"
    )

    if tmpl == TEMPLATE_SIMPLE:
        assert "#$id {color: red;background-color: blue;}" == result
    elif tmpl == TEMPLATE_SIMPLE_2:
        assert "#$id {color: red;}#$id div {background-color: blue;}" == result
