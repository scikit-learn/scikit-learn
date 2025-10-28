# Copyright (c) Povilas Kanapickas, 2019
# See LICENSE for details.

from textwrap import dedent

from twisted.trial.unittest import TestCase

from .._builder import parse_newfragment_basename, render_fragments


class TestParseNewsfragmentBasename(TestCase):
    def test_simple(self):
        """<number>.<category> generates a counter value of 0."""
        self.assertEqual(
            parse_newfragment_basename("123.feature", ["feature"]),
            ("123", "feature", 0),
        )

    def test_invalid_category(self):
        """Files without a valid category are rejected."""
        self.assertEqual(
            parse_newfragment_basename("README.ext", ["feature"]),
            (None, None, None),
        )
        self.assertEqual(
            parse_newfragment_basename("README", ["feature"]),
            (None, None, None),
        )

    def test_counter(self):
        """<number>.<category>.<counter> generates a custom counter value."""
        self.assertEqual(
            parse_newfragment_basename("123.feature.1", ["feature"]),
            ("123", "feature", 1),
        )

    def test_counter_with_extension(self):
        """File extensions are ignored."""
        self.assertEqual(
            parse_newfragment_basename("123.feature.1.ext", ["feature"]),
            ("123", "feature", 1),
        )

    def test_ignores_extension(self):
        """File extensions are ignored."""
        self.assertEqual(
            parse_newfragment_basename("123.feature.ext", ["feature"]),
            ("123", "feature", 0),
        )

    def test_non_numeric_issue(self):
        """Non-numeric issue identifiers are preserved verbatim."""
        self.assertEqual(
            parse_newfragment_basename("baz.feature", ["feature"]),
            ("baz", "feature", 0),
        )

    def test_non_numeric_issue_with_extension(self):
        """File extensions are ignored."""
        self.assertEqual(
            parse_newfragment_basename("baz.feature.ext", ["feature"]),
            ("baz", "feature", 0),
        )

    def test_dots_in_issue_name(self):
        """Non-numeric issue identifiers are preserved verbatim."""
        self.assertEqual(
            parse_newfragment_basename("baz.1.2.feature", ["feature"]),
            ("baz.1.2", "feature", 0),
        )

    def test_dots_in_issue_name_invalid_category(self):
        """Files without a valid category are rejected."""
        self.assertEqual(
            parse_newfragment_basename("baz.1.2.notfeature", ["feature"]),
            (None, None, None),
        )

    def test_dots_in_issue_name_and_counter(self):
        """Non-numeric issue identifiers are preserved verbatim."""
        self.assertEqual(
            parse_newfragment_basename("baz.1.2.feature.3", ["feature"]),
            ("baz.1.2", "feature", 3),
        )

    def test_strip(self):
        """Leading spaces and subsequent leading zeros are stripped
        when parsing newsfragment names into issue numbers etc.
        """
        self.assertEqual(
            parse_newfragment_basename("  007.feature", ["feature"]),
            ("7", "feature", 0),
        )

    def test_strip_with_counter(self):
        """Leading spaces and subsequent leading zeros are stripped
        when parsing newsfragment names into issue numbers etc.
        """
        self.assertEqual(
            parse_newfragment_basename("  007.feature.3", ["feature"]),
            ("7", "feature", 3),
        )

    def test_orphan(self):
        """Orphaned snippets must remain the orphan marker in the issue
        identifier."""
        self.assertEqual(
            parse_newfragment_basename("+orphan.feature", ["feature"]),
            ("+orphan", "feature", 0),
        )

    def test_orphan_with_number(self):
        """Orphaned snippets can contain numbers in the identifier."""
        self.assertEqual(
            parse_newfragment_basename("+123_orphan.feature", ["feature"]),
            ("+123_orphan", "feature", 0),
        )
        self.assertEqual(
            parse_newfragment_basename("+orphan_123.feature", ["feature"]),
            ("+orphan_123", "feature", 0),
        )

    def test_orphan_with_dotted_number(self):
        """Orphaned snippets can contain numbers with dots in the
        identifier."""
        self.assertEqual(
            parse_newfragment_basename("+12.3_orphan.feature", ["feature"]),
            ("+12.3_orphan", "feature", 0),
        )
        self.assertEqual(
            parse_newfragment_basename("+orphan_12.3.feature", ["feature"]),
            ("+orphan_12.3", "feature", 0),
        )

    def test_orphan_all_digits(self):
        """Orphaned snippets can consist of only digits."""
        self.assertEqual(
            parse_newfragment_basename("+123.feature", ["feature"]),
            ("+123", "feature", 0),
        )


class TestNewsFragmentsOrdering(TestCase):
    """
    Tests to ensure that issues are ordered correctly in the output.

    This tests both ordering of issues within a fragment and ordering of
    fragments within a section.
    """

    template = dedent(
        """
    {% for section_name, category in sections.items() %}
    {% if section_name %}# {{ section_name }}{% endif %}
    {%- for category_name, issues in category.items() %}
    ## {{ category_name }}
    {% for issue, numbers in issues.items() %}
    - {{ issue }}{% if numbers %} ({{ numbers|join(', ') }}){% endif %}

    {% endfor %}
    {% endfor -%}
    {% endfor -%}
    """
    )

    def render(self, fragments):
        return render_fragments(
            template=self.template,
            issue_format=None,
            fragments=fragments,
            definitions={},
            underlines=[],
            wrap=False,
            versiondata={},
        )

    def test_ordering(self):
        """
        Issues are ordered first by the non-text component, then by their number.

        For backwards compatibility, issues with no number are grouped first and issues
        which are only a number are grouped last.

        Orphan news fragments are always last, sorted by their text.
        """
        output = self.render(
            {
                "": {
                    "feature": {
                        "Added Cheese": ["10", "gh-25", "gh-3", "4"],
                        "Added Fish": [],
                        "Added Bread": [],
                        "Added Milk": ["gh-1"],
                        "Added Eggs": ["gh-2", "random"],
                    }
                }
            },
        )
        # "Eggs" are first because they have an issue with no number, and the first
        # issue for each fragment is what is used for sorting the overall list.
        assert output == dedent(
            """
            ## feature
            - Added Eggs (random, gh-2)
            - Added Milk (gh-1)
            - Added Cheese (gh-3, gh-25, #4, #10)
            - Added Bread
            - Added Fish
"""
        )
