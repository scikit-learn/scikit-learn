from contextlib import suppress

from docutils import nodes
from docutils.parsers.rst import Directive

from sklearn.utils import all_estimators
from sklearn.utils._test_common.instance_generator import _construct_instances
from sklearn.utils._testing import SkipTest


class AllowNanEstimators(Directive):
    @staticmethod
    def make_paragraph_for_estimator_type(estimator_type):
        intro = nodes.list_item()
        intro += nodes.strong(text="Estimators that allow NaN values for type ")
        intro += nodes.literal(text=f"{estimator_type}")
        intro += nodes.strong(text=":\n")
        exists = False
        lst = nodes.bullet_list()
        for name, est_class in all_estimators(type_filter=estimator_type):
            with suppress(SkipTest):
                # Here we generate the text only for one instance. This directive
                # should not be used for meta-estimators where tags depend on the
                # sub-estimator.
                est = next(_construct_instances(est_class))

                if est.__sklearn_tags__().input_tags.allow_nan:
                    module_name = ".".join(est_class.__module__.split(".")[:2])
                    class_title = f"{est_class.__name__}"
                    class_url = f"./generated/{module_name}.{class_title}.html"
                    item = nodes.list_item()
                    para = nodes.paragraph()
                    para += nodes.reference(
                        class_title, text=class_title, internal=False, refuri=class_url
                    )
                    exists = True
                    item += para
                    lst += item
        intro += lst
        return [intro] if exists else None

    def run(self):
        lst = nodes.bullet_list()
        for i in ["cluster", "regressor", "classifier", "transformer"]:
            item = self.make_paragraph_for_estimator_type(i)
            if item is not None:
                lst += item
        return [lst]


def setup(app):
    app.add_directive("allow_nan_estimators", AllowNanEstimators)

    return {
        "version": "0.1",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
