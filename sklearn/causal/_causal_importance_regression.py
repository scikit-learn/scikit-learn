from ._causal_importance import CausalFeatureImportance

class CausalFeatureImportanceRegression(CausalFeatureImportance):
    def __init__(self,model,treatment,control,label):
        super().__init__(model,treatment,control,label)
        
