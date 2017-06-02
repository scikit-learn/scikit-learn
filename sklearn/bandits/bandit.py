# -*- coding: utf-8 -*-

"""
"""
from abc import abstractmethod, ABCMeta

# Author: Eustache Diemert <eustache@diemert.fr>
# License: BSD Style.

class Bandit(object):
    """Abstract class common to all bandits"""

    __metaclass__ = ABCMeta

    @abstractmethod
    def recommend(self, visitor, possible_actions):
        """Chooses an action among possible_actions.

        Parameters
        ----------
        visitor : hashable, identify visitors to which we recommend an
                  action (may not be used at all by the implementation)
        possible_actions : array-like, identify all possible actions
                           we can recommend to the visitor.

        Returns
        -------
        action : one of possible_actions or None
        """

    @abstractmethod
    def update(self, visitor, action, reward):
        """Gives feedback to the algorithm about the result of the
        action we recommended to a visitor.

        Parameters
        ----------
        visitor : hashable, identify visitor
        action : action the algorithm recommended via recommend()
        reward : binary (0/1) or float [0,1] reward. 1 means best reward.

        Returns
        -------
        None
        """
