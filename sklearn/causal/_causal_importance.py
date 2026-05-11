import numpy as np
from utils import clone_model

class CausalFeatureImportance:
    def __init__(self,model,treatment,control,label):
        self.model=model
        self.treatment=treatment
        self.control=control
        self.label=label
