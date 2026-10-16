# SVM Sentiment Analysis Example

This example demonstrates three-class sentiment analysis using
scikit-learn.

The classifier predicts whether a sentence is:

- positive
- negative
- neutral

The implementation uses a pipeline containing:

1. `TfidfVectorizer` for converting text into numerical features
2. `LinearSVC` for Support Vector Machine classification

## Run the example

From the root of the scikit-learn repository:

```bash
python examples/text/sentiment_analysis_svm.py
