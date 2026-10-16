# scikit-learn BaseEstimator fix for feature names
class BaseEstimator:
    def to(self, device):
        """
        Move estimator to specified device.
        Preserves feature_names_in_ attribute.
        """
        # Move main estimator data
        moved = self.__class__(**self.get_params())
        
        # Fix: preserve feature_names_in_ when moving estimator
        if hasattr(self, 'feature_names_in_'):
            moved.feature_names_in_ = self.feature_names_in_
        
        return moved
