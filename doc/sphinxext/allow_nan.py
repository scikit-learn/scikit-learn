from sklearn.utils import all_estimators
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import OrdinalEncoder, OneHotEncoder
from sklearn.compose import ColumnTransformer
from io import StringIO
from docutils import nodes

from docutils.parsers.rst import Directive


class Allow_Nan(Directive):
    def run(self):
        names1 = []
        names2 = []
        names3 = []
        names4 = []
        output = StringIO()
        output.write(
            f"List of estimators that allow NaN values for type_filter=cluster : "
        )
        for name, Estimator in all_estimators(type_filter="cluster"):
            if Estimator()._get_tags().get("allow_nan") == True:
                names1.append(name)
        output.write(f" {names1} ")

        output2 = StringIO()
        output2.write(
            f" List of estimators that allow NaN values for type_filter=regressor : "
        )
        for name, Estimator in all_estimators(type_filter="regressor"):
            try:
                if Estimator()._get_tags().get("allow_nan") == True:
                    names2.append(name)
            except TypeError:
                base = LogisticRegression(solver="lbfgs")
                if Estimator(base)._get_tags().get("allow_nan") == True:
                    names2.append(name)
        output2.write(f" {names2} ")

        output3 = StringIO()
        output3.write(
            f" List of estimators that allow NaN values for type_filter=classifier : "
        )
        for name, Estimator in all_estimators(type_filter="classifier"):
            try:
                if Estimator()._get_tags().get("allow_nan") == True:
                    names3.append(name)
            except TypeError:
                base = LogisticRegression(solver="lbfgs")
                if Estimator(base)._get_tags().get("allow_nan") == True:
                    names3.append(name)
        output3.write(f" {names3} ")

        output4 = StringIO()
        output4.write(
            f" List of estimators that allow NaN values for type_filter=transformer : "
        )
        for name, Estimator in all_estimators(type_filter="transformer"):
            try:
                if Estimator()._get_tags().get("allow_nan") == True:
                    names4.append(name)
            except TypeError:
                if name == "ColumnTransformer":
                    if (
                        ColumnTransformer(
                            [
                                ("ordinal", OrdinalEncoder(), [0, 1]),
                                ("nominal", OneHotEncoder(), [2, 3]),
                            ]
                        )
                        ._get_tags()
                        .get("allow_nan")
                        == True
                    ):
                        names4.append(name)
        output4.write(f" {names4} ")

        output = output.getvalue()
        output2 = output2.getvalue()
        output3 = output3.getvalue()
        output4 = output4.getvalue()
        paragraph_node = nodes.paragraph(text=output)
        paragraph_node2 = nodes.paragraph(text=output2)
        paragraph_node3 = nodes.paragraph(text=output3)
        paragraph_node4 = nodes.paragraph(text=output4)

        return [paragraph_node, paragraph_node2, paragraph_node3, paragraph_node4]


def setup(app):

    app.add_directive("allow_nan", Allow_Nan)

    return {
        "version": "0.1",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
