import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._visible import VisibleValidator
    from ._valuessrc import ValuessrcValidator
    from ._values import ValuesValidator
    from ._tickvalssrc import TickvalssrcValidator
    from ._tickvals import TickvalsValidator
    from ._ticktextsrc import TicktextsrcValidator
    from ._ticktext import TicktextValidator
    from ._tickformat import TickformatValidator
    from ._templateitemname import TemplateitemnameValidator
    from ._range import RangeValidator
    from ._name import NameValidator
    from ._multiselect import MultiselectValidator
    from ._label import LabelValidator
    from ._constraintrange import ConstraintrangeValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._visible.VisibleValidator",
            "._valuessrc.ValuessrcValidator",
            "._values.ValuesValidator",
            "._tickvalssrc.TickvalssrcValidator",
            "._tickvals.TickvalsValidator",
            "._ticktextsrc.TicktextsrcValidator",
            "._ticktext.TicktextValidator",
            "._tickformat.TickformatValidator",
            "._templateitemname.TemplateitemnameValidator",
            "._range.RangeValidator",
            "._name.NameValidator",
            "._multiselect.MultiselectValidator",
            "._label.LabelValidator",
            "._constraintrange.ConstraintrangeValidator",
        ],
    )
