"""
====================================================================
Custom Tokenization for Non-European Languages (e.g. Hindi, Marathi)
====================================================================

The default tokenizer in :class:`~sklearn.feature_extraction.text.CountVectorizer`
is broadly engineered for languages with scripts where letters are typed sequentially
and don't rely heavily on separate floating modifiers. It uses the regular expression
`r"(?u)\\b\\w\\w+\\b"`.

In Python's `re` module, `\\w` matches letters and numbers, but it does not match
with a Unicode nonspacing mark (Category Mn). For some scripts like Devanagari,
the vowel modifiers (maatras) are part of this category. Consequently, the default
tokenizer treats these modifiers as punctuation or delimiters splitting single words
into meaningless fragments.

This example ilustrates how to override the default `token_pattern` to correctly
tokenize text in scripts, that include nonspacing marks like Devanagari, by
explicitly including its Unicode block.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Basic Case: An atomic example
# -----------------------------
# Using a single word with a modifier. The word "हिंदी" (Hindi) contains modifiers
# that are not captured by the standard `\w` character class.

from sklearn.feature_extraction.text import CountVectorizer

basic_text = "हिंदी"

# We are using `build_analyzer()` to view the tokens directly since the default
# analyzer completely destroys the word, filtering it out entirely because the
# fractured fragments are less than 2 characters long.
vectorizer_default = CountVectorizer()
analyzer_default = vectorizer_default.build_analyzer()
print("Default Tokenizer (Basic Text):", analyzer_default(basic_text))

# Passing a custom token_pattern that includes the Devanagari Unicode block
# (\u0900-\u097F), the base characters and their modifiers are captured together.
vectorizer_devanagari = CountVectorizer(token_pattern=r"(?u)[\w\u0900-\u097F]+")
analyzer_devanagari = vectorizer_devanagari.build_analyzer()
print("Custom Tokenizer (Basic Text): ", analyzer_devanagari(basic_text))

# %%
# Complex Case: A complete sentence
# ---------------------------------
# Using the same sentence user rachanagusain stated in the issue #20211 it contains
# modifiers, numbers and punctuations. The default tokenizer strips out the modifiers
# and leaves broken fragments of words.

complex_text = [
    "आई- जाइयै नुहाड़ा ध्यान उस्सै भेठा जा’रदा हा।",  # noqa: RUF001
    "य बात सुणते ही सब गौं वाल नजदीक जंगलै तरफ भाज् और बांनर पकड़ पकड़ बेर 100- 100 रु में बेचंण लाग्।",
]

vectorizer_default.fit(complex_text)
print("\nDefault Tokenizer (Complex):\n", vectorizer_default.get_feature_names_out())


# The custom pattern for Devanagari ignores the punctuations and whitespace while
# keeping the full words and numbers intact.
vectorizer_devanagari.fit(complex_text)
print("\nCustom Tokenizer (Complex):\n", vectorizer_devanagari.get_feature_names_out())
