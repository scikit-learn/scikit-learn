"""
Three-class sentiment analysis using TF-IDF and a linear SVM.

This example predicts whether a sentence is positive, negative, or neutral.
"""

import argparse

from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics import classification_report
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.svm import LinearSVC


TEXTS = [
    "I love this product",
    "This product is excellent",
    "The service was wonderful",
    "I am very happy with the result",
    "This movie was amazing",
    "The experience was fantastic",
    "The product works perfectly",
    "I really enjoyed this service",
    "This is a great application",
    "I hate this product",
    "This product is terrible",
    "The service was awful",
    "I am very disappointed",
    "This movie was boring",
    "The experience was horrible",
    "The product stopped working",
    "I regret buying this product",
    "This service is extremely poor",
    "The package arrived today",
    "The meeting starts at noon",
    "The product is available online",
    "I received the email this morning",
    "The movie is two hours long",
    "The store opens at nine",
    "The application was released yesterday",
    "The report contains ten pages",
    "The train arrives at five",
    "The appointment is scheduled for Monday",
]

LABELS = [
    "positive",
    "positive",
    "positive",
    "positive",
    "positive",
    "positive",
    "positive",
    "positive",
    "positive",
    "negative",
    "negative",
    "negative",
    "negative",
    "negative",
    "negative",
    "negative",
    "negative",
    "negative",
    "neutral",
    "neutral",
    "neutral",
    "neutral",
    "neutral",
    "neutral",
    "neutral",
    "neutral",
    "neutral",
    "neutral",
]


def build_model():
    """Create a TF-IDF and linear SVM classification pipeline."""
    return make_pipeline(
        TfidfVectorizer(
            lowercase=True,
            ngram_range=(1, 2),
        ),
        LinearSVC(),
    )


def train_model():
    """Train the sentiment classifier and return the fitted model."""
    model = build_model()
    model.fit(TEXTS, LABELS)
    return model


def main():
    """Train the model, evaluate it, and optionally classify new text."""
    parser = argparse.ArgumentParser(
        description="Classify text as positive, negative, or neutral."
    )
    parser.add_argument(
        "--text",
        help="A sentence to classify after training.",
    )
    args = parser.parse_args()

    x_train, x_test, y_train, y_test = train_test_split(
        TEXTS,
        LABELS,
        test_size=0.3,
        random_state=42,
        stratify=LABELS,
    )

    model = build_model()
    model.fit(x_train, y_train)

    predictions = model.predict(x_test)

    print("Model evaluation:")
    print(
        classification_report(
            y_test,
            predictions,
            labels=["negative", "neutral", "positive"],
            zero_division=0,
        )
    )

    if args.text:
        final_model = train_model()
        prediction = final_model.predict([args.text])[0]
        print(f"Text: {args.text}")
        print(f"Predicted sentiment: {prediction}")


if __name__ == "__main__":
    main()
