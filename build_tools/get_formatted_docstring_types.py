"""Helper script to get the docstring type of estimator or function."""
from importlib import import_module
import argparse
import inspect


from sklearn.utils._typing import get_docstring_annotations


parser = argparse.ArgumentParser(
    description=("Generates typed docstring for a specific scikit-learn "
                 "class or function"))
parser.add_argument('object', help=("scikit-learn object, for example "
                                    "linear_model.LogisticRegression"))

args = parser.parse_args()
object_input = args.object
object_split = object_input.split(".")

module = "sklearn." + ".".join(object_split[:-1])
obj_str = object_split[-1]
obj = getattr(import_module(module), obj_str)

print("Parameters")
print("----------")
if inspect.isclass(obj):
    formatted_annotations = get_docstring_annotations(obj.__init__)
else:  # function
    formatted_annotations = get_docstring_annotations(obj)
for name, annotation in formatted_annotations.items():
    print(f"{name} : {annotation}")
