# Author: Arturo Amor
# License: 3-clause BSD
"""Test the example recommender plugin."""

import codecs
import os
import re

import pytest

import sphinx_gallery.gen_rst as sg
from sphinx_gallery.recommender import ExampleRecommender, _write_recommendations


def test_example_recommender_methods():
    """Test dict_vectorizer and compute_tf_idf methods."""
    np = pytest.importorskip("numpy")
    recommender = ExampleRecommender()
    D = [{"foo": 1, "bar": 2}, {"foo": 3, "baz": 1}]

    X = recommender.dict_vectorizer(D)
    expected_X = np.array([[2.0, 0.0, 1.0], [0.0, 1.0, 3.0]])
    np.testing.assert_array_equal(X, expected_X)

    X_tfidf = recommender.compute_tf_idf(X)
    expected_X_tfidf = np.array(
        [[0.94215562, 0.0, 0.33517574], [0.0, 0.42423963, 0.90554997]]
    )
    np.testing.assert_array_almost_equal(X_tfidf, expected_X_tfidf)


def test_recommendation_files(gallery_conf):
    """Test generated files and that recommendations are relevant."""
    pytest.importorskip("numpy")
    gallery_conf["recommender"].update(
        [("enable", True), ("rubric_header", "Custom header")]
    )
    file_dict = {
        "fox_jumps_dog.py": "The quick brown fox jumped over the lazy dog",
        "dog_sleeps.py": "The lazy dog slept all day",
        "fox_eats_dog_food.py": "The quick brown fox ate the lazy dog's food",
        "dog_jumps_fox.py": "The quick dog jumped over the lazy fox",
    }

    for file_name, content in file_dict.items():
        file_path = os.path.join(gallery_conf["gallery_dir"], file_name)
        with open(file_path, "w") as f:
            f.write(content)
        sg.save_rst_example("example_rst", file_path, 0, 0, gallery_conf)

        test_file = re.sub(r"\.py$", ".rst", file_path)
        recommendation_file = re.sub(r"\.py$", ".recommendations", file_name)
        with codecs.open(test_file) as f:
            rst = f.read()

        assert recommendation_file in rst

    py_files = [
        fname
        for fname in os.listdir(gallery_conf["gallery_dir"])
        if os.path.splitext(fname)[1] == ".py"
    ]
    gallery_py_files = [
        os.path.join(gallery_conf["gallery_dir"], fname) for fname in py_files
    ]
    recommender = ExampleRecommender(n_examples=1, min_df=1)
    recommender.fit(gallery_py_files)
    recommended_example = recommender.predict(file_path)  # dog_jumps_fox.py

    assert os.path.basename(recommended_example[0]) == "fox_jumps_dog.py"

    # _write_recommendations needs a thumbnail, for writing the
    # `_thumbnail_div` we then create a blank png
    thumb_path = os.path.join(gallery_conf["gallery_dir"], "images/thumb")
    os.makedirs(thumb_path, exist_ok=True)
    png_file = "sphx_glr_fox_jumps_dog_thumb.png"
    png_file_path = os.path.join(thumb_path, png_file)
    with open(png_file_path, "wb") as f:
        b"\x89PNG\r\n\x1a\n"  # generic png file signature

    recommendation_file = re.sub(r"\.py$", ".recommendations", file_path)
    _write_recommendations(recommender, file_path, gallery_conf)
    with codecs.open(recommendation_file) as f:
        rst = f.read()
    assert ".. rubric:: Custom header" in rst
