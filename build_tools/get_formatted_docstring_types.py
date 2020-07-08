"""Helper script to get the docstring type of estimator or function."""
from importlib import import_module
import argparse
import inspect


from sklearn.utils._typing import get_annotations


parser = argparse.ArgumentParser(
    description=("Generates typed docstring for a specific scikit-learn "
                 "class or function"))
parser.add_argument('object', help=("scikit-learn object, for example "
                                    "linear_model.LogisticRegression"))

args = parser.parse_args()
object_input = args.object
object_split = object_input.split(".")

module = "sklearn." + ".".join(object_split[:-1])
instance_str = object_split[-1]
instance = getattr(import_module(module), instance_str)


print("Parameters")
print("----------")
if inspect.isclass(instance):
    formatted_annotations = get_annotations(instance.__init__)
else:
    formatted_annotations = get_annotations(instance)
for name, annotation in formatted_annotations.items():
    print(f"{name} : {annotation}")


if inspect.isclass(instance):
    print()
    print("Attributes")
    print("----------")
    formatted_annotations = get_annotations(instance)
    for name, annotation in formatted_annotations.items():
        print(f"{name} : {annotation}")
