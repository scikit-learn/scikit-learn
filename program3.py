
def clean_text(text):
    """Normalize text before passing it to an NLP model."""
    return " ".join(text.lower().strip().split())


if __name__ == "__main__":
    sample = "  Machine   Learning PROJECT  "
    print(clean_text(sample))
