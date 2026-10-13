from examples.text.sentiment_analysis_svm import (
    LABELS,
    TEXTS,
    build_model,
    train_model,
)


def test_texts_and_labels_have_same_length():
    assert len(TEXTS) == len(LABELS)


def test_model_can_predict_three_sentiments():
    model = train_model()

    predictions = model.predict(
        [
            "This product is excellent",
            "This service is terrible",
            "The meeting starts at noon",
        ]
    )

    assert predictions[0] == "positive"
    assert predictions[1] == "negative"
    assert predictions[2] == "neutral"


def test_model_pipeline_contains_tfidf_and_svm():
    model = build_model()

    assert "tfidfvectorizer" in str(model).lower()
    assert "linearsvc" in str(model).lower()
