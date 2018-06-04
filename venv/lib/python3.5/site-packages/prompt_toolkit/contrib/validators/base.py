from __future__ import unicode_literals
from prompt_toolkit.validation import Validator, ValidationError
from six import string_types


class SentenceValidator(Validator):
    """
    Validate input only when it appears in this list of sentences.

    :param sentences: List of sentences.
    :param ignore_case: If True, case-insensitive comparisons.
    """
    def __init__(self, sentences, ignore_case=False, error_message='Invalid input', move_cursor_to_end=False):
        assert all(isinstance(s, string_types) for s in sentences)
        assert isinstance(ignore_case, bool)
        assert isinstance(error_message, string_types)

        self.sentences = list(sentences)
        self.ignore_case = ignore_case
        self.error_message = error_message
        self.move_cursor_to_end = move_cursor_to_end

        if ignore_case:
            self.sentences = set([s.lower() for s in self.sentences])

    def validate(self, document):
        if document.text not in self.sentences:
            if self.move_cursor_to_end:
                index = len(document.text)
            else:
                index = 0

            raise ValidationError(cursor_position=index,
                                  message=self.error_message)
