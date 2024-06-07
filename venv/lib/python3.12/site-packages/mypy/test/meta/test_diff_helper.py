import io

from mypy.test.helpers import Suite, diff_ranges, render_diff_range


class DiffHelperSuite(Suite):
    def test_render_diff_range(self) -> None:
        expected = ["hello", "world"]
        actual = ["goodbye", "world"]

        expected_ranges, actual_ranges = diff_ranges(expected, actual)

        output = io.StringIO()
        render_diff_range(expected_ranges, expected, output=output)
        assert output.getvalue() == "  hello (diff)\n  world\n"
        output = io.StringIO()
        render_diff_range(actual_ranges, actual, output=output)
        assert output.getvalue() == "  goodbye (diff)\n  world\n"

        expected = ["a", "b", "c", "d", "e", "f", "g", "h", "circle", "i", "j"]
        actual = ["a", "b", "c", "d", "e", "f", "g", "h", "square", "i", "j"]

        expected_ranges, actual_ranges = diff_ranges(expected, actual)

        output = io.StringIO()
        render_diff_range(expected_ranges, expected, output=output, indent=0)
        assert output.getvalue() == "a\nb\nc\n...\nf\ng\nh\ncircle (diff)\ni\nj\n"
        output = io.StringIO()
        render_diff_range(actual_ranges, actual, output=output, indent=0)
        assert output.getvalue() == "a\nb\nc\n...\nf\ng\nh\nsquare (diff)\ni\nj\n"

    def test_diff_ranges(self) -> None:
        a = ["hello", "world"]
        b = ["hello", "world"]

        assert diff_ranges(a, b) == (
            [(0, 0), (0, 2), (2, 2), (2, 2)],
            [(0, 0), (0, 2), (2, 2), (2, 2)],
        )

        a = ["hello", "world"]
        b = ["goodbye", "world"]

        assert diff_ranges(a, b) == (
            [(0, 1), (1, 2), (2, 2), (2, 2)],
            [(0, 1), (1, 2), (2, 2), (2, 2)],
        )
