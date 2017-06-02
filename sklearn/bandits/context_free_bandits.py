# -*- coding: utf-8 -*-

"""
The :mod: sklearn.bandits.context_free_bandits module implements
Context-Free Multi-Armed Bandits. These are online learning methods related to
recommender systems and reinforcement learning.
See also http://en.wikipedia.org/wiki/Multi-armed_bandit
"""

# Author: Eustache Diemert <eustache@diemert.fr>
# License: BSD Style.

import random

class Bandit(object):
        """Abstract class common to all bandits"""

        @abstractmethod
        def recommend(self, visitor, possible_actions):
            """Chooses an action among possible_actions.

            Parameters
            ----------
            visitor : hashable, identify visitors to which we recommend an
                      action (may not be used at all y the implementation)
            possible_actions : array-like, identify all possible actions
                               we can recommend to the visitor.

            Returns
            -------
            action : one of possible_actions
            """

        @abstractmethod
        def update(self, visitor, action, reward):
            """Gives feedback to the algorithm about the result of the
            action we recommended to visitor.

            Parameters
            ----------
            visitor : hashable, identify visitor
            action : action the algorithm recommended via recommend()
            reward : binary (0/1) or float [0,1] reward. 1 means best reward.

            Returns
            -------
            None
            """

class EpsilonGreedyBandit(Bandit):
    """
    The best action (as much as the algorithm knows so far) is selected for
    a proportion 1 - \epsilon of the trials, and another action is randomly
    selected (with uniform probability) for a proportion \epsilon.

    Parameters
    ----------
    epsilon : float in [0,1], proportion of trials in which we select a
              random action
    
    """

    def __init__(self, epsilon=0.1):
        # raises TypeError if conversion not supported
        self.epsilon = float(epsilon)
        if self.epsilon < 0.0 or epsilon > 1.0:
            raise ValueError('epsilon must be in range [0,1]')
        self.actions_played = dict()
        self.actions_rewarded = dict()

    def _expected_reward(self, action):
        if action not in self.actions_played:
            raise ValueError('action %s never played'%repr(action))
        mean_reward = self.actions_rewarded.get(action,0.0)
        mean_reward /= self.actions_played[action]
        return mean_reward

    def recommend(self, visitor, possible_actions):
        if random.uniform(0.0,1.0) < self.epsilon:
            return random.choice(possible_actions)
        # sort actions by expected reward and selects best
        # mean rewarded action
        action_rewards = [ (self._expected_reward(action),action) for
                           action in self.actions_played.iterkeys() ]
        action_rewards.sort(reverse=True)
        chosen_action = action_rewards[0][1]
        self.actions_played[chosen_action] = self.actions_played.get(chosen_action,
                                                                     0.0) + 1.0
        return chosen_action

    def update(self, visitor, action, reward):
        # raises TypeError if conversion not supported
        reward = float(reward)
        if reward < 0.0 or reward > 1.0:
            raise ValueError('reward must be in range [0,1]')
        self.actions_rewarded[action] = self.actions_rewarded.get(action,0.0) + 1.0
