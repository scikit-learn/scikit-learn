import runpy
import types
import numpy as np
import sys
from pathlib import Path


def _load_export_text():
    root = Path(__file__).resolve().parent.parent
    skl_path = root / 'sklearn'

    sklearn_pkg = types.ModuleType('sklearn')
    sklearn_pkg.__path__ = [str(skl_path)]
    sys_mod = {'sklearn': sklearn_pkg}

    utils_pkg = types.ModuleType('sklearn.utils')
    utils_pkg.__path__ = [str(skl_path / 'utils')]
    sys_mod['sklearn.utils'] = utils_pkg

    validation_mod = types.ModuleType('sklearn.utils.validation')
    validation_mod.check_is_fitted = lambda est, attr: None
    sys_mod['sklearn.utils.validation'] = validation_mod

    tree_pkg = types.ModuleType('sklearn.tree')
    tree_pkg.__path__ = [str(skl_path / 'tree')]
    tree_pkg.DecisionTreeClassifier = type('DecisionTreeClassifier', (), {})
    tree_pkg.DecisionTreeRegressor = type('DecisionTreeRegressor', (), {})
    sys_mod['sklearn.tree'] = tree_pkg

    criterion_mod = types.ModuleType('sklearn.tree._criterion')
    sys_mod['sklearn.tree._criterion'] = criterion_mod

    tree_mod = types.ModuleType('sklearn.tree._tree')
    tree_mod.TREE_UNDEFINED = -2
    sys_mod['sklearn.tree._tree'] = tree_mod

    for name, mod in sys_mod.items():
        sys.modules[name] = mod

    module = runpy.run_module('sklearn.tree.export', run_name='__main__')
    return module['export_text']


def test_export_text_single_feature():
    export_text = _load_export_text()

    class Tree:
        n_features = 1
        feature = np.array([0, -2, -2])
        threshold = np.array([0.5, -2, -2])
        children_left = np.array([1, -1, -1])
        children_right = np.array([2, -1, -1])
        value = np.array([[[1, 1]], [[1, 0]], [[0, 1]]], dtype=float)
        n_outputs = 1
        n_classes = [2]

    class Estimator:
        tree_ = Tree()
        classes_ = np.array([0, 1])

    report = export_text(Estimator(), feature_names=['a'])
    assert report.startswith('|--- a')