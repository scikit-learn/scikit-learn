from sklearn.utils import all_estimators
from sklearn.compose import ColumnTransformer
from sklearn.utils.estimator_checks import _construct_instance
from io import StringIO
from docutils import nodes
import warnings

from docutils.parsers.rst import Directive


class Allow_Nan(Directive):
    @staticmethod
    def make_paragraph_for_estimator_type(estimator_type):
        output = StringIO()
        output.write(
            f"* List of estimators that allow NaN values for type *{estimator_type}*:\n"
        )

        exists = False
        for name, est_class in all_estimators(type_filter=estimator_type):
            try:
                est = _construct_instance(est_class)
            except:
                if name == "ColumnTransformer":
                    est = ColumnTransformer(None)
                else:
                    warnings.warn(
                        f"Estimator {est_class.__name__} failed to construct."
                    )

            if est._get_tags().get("allow_nan"):
                module_name = est_class.__module__
                class_name = est_class.__name__
                output.write(f" * :class:`{module_name}.{class_name}`\n")
                exists = True
        return nodes.paragraph(text=output.getvalue()) if exists else None

    def run(self):
        output = []
        for i in ["cluster", "regressor", "classifier", "transformer"]:
            paragraph = self.make_paragraph_for_estimator_type(i)
            if paragraph is not None:
                output.append(paragraph)
        return output


def setup(app):

    app.add_directive("allow_nan", Allow_Nan)

    return {
        "version": "0.1",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
