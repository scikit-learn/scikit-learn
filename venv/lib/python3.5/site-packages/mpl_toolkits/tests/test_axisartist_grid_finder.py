from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from mpl_toolkits.axisartist.grid_finder import (
    FormatterPrettyPrint,
    MaxNLocator)


def test_pretty_print_format():
    locator = MaxNLocator()
    locs, nloc, factor = locator(0, 100)

    fmt = FormatterPrettyPrint()

    assert fmt("left", None, locs) == \
        [r'$\mathdefault{%d}$' % (l, ) for l in locs]
