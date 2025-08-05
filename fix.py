# fixes/fix_scikit-learn_26711.py

class LabelEncoder:
    def set_output(self, transform=None):
        pass  # placeholder for compatibility
