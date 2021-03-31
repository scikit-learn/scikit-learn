"""
The :mod:`sklearn.domain_adaptation` module implements domain adaptation
algorithms. These algorithms minimize the domain gaps between domains
to improve the classification performance. This module includes Transfer 
Component Analysis (TCA) and Balanced Distribution Adaptation (BDA).
"""

from ._transfer_component_analysis import TransferComponentAnalysis
from ._balanced_distribution_adaptation import BalancedDistributionAdaptation

__all__ = ['TransferComponentAnalysis', 'BalancedDistributionAdaptation']