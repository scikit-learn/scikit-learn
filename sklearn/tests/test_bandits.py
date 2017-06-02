# -*- coding: utf-8 -*-

# Author: Eustache Diemert <eustache@diemert.fr>
# License: BSD Style.

import math
import nose
from nose.tools import timed
import numpy as np
from numpy.oldnumeric.random_array import binomial

from sklearn.bandits.context_free_bandits import EpsilonGreedyBandit

#############################################################################
# Helper Classes to generate datasets

# set random seed to get reproduceable tests
np.random.seed(0)

class RandomStaticExperiment(object):
    '''
    Generates a random experiment which can be used for testing purposes. 
    Each possible action is drawn from a binomial with reward probability p 
    '''

    def __init__(self, duration, reward_probabilities):
        self.duration = int(duration)
        self.expected_reward = {}
        self.mean_reward = {}
        self.data = {}
        arm = 0
        for p in reward_probabilities:
            self.data[arm] = binomial(1, p, self.duration)
            self.expected_reward[arm] = p
            self.mean_reward[arm] = sum(self.data[arm]) / float(self.duration)
            print "generated arm: %s with expect. reward: %.6f (real: %.6f)" % (arm, p, self.mean_reward[arm])
            arm += 1

    def __iter__(self):
        self.choices = {}
        return self

    def next(self):
        self.duration -= 1
        if self.duration >= 0:
            return [ (arm, self.data[arm][self.duration]) for arm in self.data.iterkeys() ]
        raise StopIteration

    def __call__(self, recommender):
        visitor = 0
        rewards_sum = 0
        for time_step in self:
            actions_rewards = dict(time_step)
            possible_actions, _rewards = zip(*time_step)
            chosen_action = recommender.recommend(visitor, possible_actions)
            reward = actions_rewards[chosen_action]
            rewards_sum += reward
            recommender.update(visitor, chosen_action, reward)
            self.choices[chosen_action] = self.choices.get(chosen_action, 0) + 1
        print "choices made:", self.choices
        return self.evaluate(rewards_sum)

    def evaluate(self, realized_rewards_sum):
        #best_action = None
        best_reward_sum = -1.0
        for _action, rewards in self.data.iteritems():
            rewards_sum = sum(rewards)
            if best_reward_sum < rewards_sum:
                #best_action = _action
                best_reward_sum = rewards_sum
        print "evaluation: realized %d best %d (%.6f %%)" % (realized_rewards_sum,
                                                             best_reward_sum,
                                                             100 * realized_rewards_sum / float(best_reward_sum))
        return (realized_rewards_sum, best_reward_sum)

class RandomDynamicExperiment(object):
    '''
    Generates a random experiment which can be used for testing bandits robustness 
    against a dynamic adversary. 
    Each possible action is drawn from a binomial with reward probability p(t)
    where t is a timestep (in other words, the rewards probabilities change over time) 
    '''
    
    def __init__(self, step_duration, reward_probabilities):
        """
        @param step_duration : int
        @param reward_probabilities : numpy ndarray of shape [k,s] representing 
               the K arms and their mean reward during the S steps
        """
        self.step_duration = step_duration
        self.reward_probabilities = reward_probabilities
        self.mean_reward = {}
        self.nb_arms = reward_probabilities.shape[0]
        self.nb_steps = reward_probabilities.shape[1]
        self.data = np.ndarray(shape=(self.nb_arms,
                                      self.nb_steps * step_duration),
                                      dtype=int)
        for arm in range(0, self.nb_arms):
            for step in range(0, self.nb_steps):
                d = binomial(1, reward_probabilities[arm, step], self.step_duration)
                self.data[arm][step * step_duration:(step + 1) * step_duration] = d
            self.mean_reward[arm] = sum(self.data[arm]) / float(self.step_duration)
            print "generated arm: %s with expt. rewards %s X %d (real reward: %.6f) " % (arm,
                                                                                        reward_probabilities[arm],
                                                                                        self.step_duration,
                                                                                        self.mean_reward[arm])
            arm += 1

    def __iter__(self):
        self.current_index = -1
        self.choices = {}
        return self

    def next(self):
        self.current_index += 1
        if self.current_index >= self.nb_steps * self.step_duration:
            raise StopIteration
        return [ (arm, self.data[arm][self.current_index]) for arm in range(0, self.nb_arms) ]

    def __call__(self, recommender):
        visitor = 0
        rewards_sum = 0
        for time_step in self:
            actions_rewards = dict(time_step)
            possible_actions, _rewards = zip(*time_step)
            chosen_action = recommender.recommend(visitor, possible_actions)
            reward = actions_rewards[chosen_action]
            rewards_sum += reward
            recommender.update(visitor, chosen_action, reward)
            self.choices[chosen_action] = self.choices.get(chosen_action, 0) + 1
        print "choices made:", self.choices
        return self.evaluate(rewards_sum)

    def evaluate(self, realized_rewards_sum):
        #best_action = None
        best_reward_sum = -1.0
        for rewards in self.data:
            rewards_sum = sum(rewards)
            if best_reward_sum < rewards_sum:
                #best_action = _action
                best_reward_sum = rewards_sum
        print "evaluation: realized %d best %d (%.6f %%)" % (realized_rewards_sum,
                                                             best_reward_sum,
                                                             100 * realized_rewards_sum / float(best_reward_sum))
        return (realized_rewards_sum, best_reward_sum)

#############################################################################
# The tests


# create a static adversary with 6 arms and duration=1000 iterations
basic_static_experiment = RandomStaticExperiment(1000, [0.2,
                                                        0.6,
                                                        0.1,
                                                        0.2,
                                                        0.3,
                                                        0.5])

# should pass in < 50ms per call
@timed(basic_static_experiment.duration*0.05)
def test_epsilon_greedy():
    recommender = EpsilonGreedyBandit()
    realized_rewards_sum, best_reward_sum = basic_static_experiment(recommender)
    realized_vs_best_ratio = math.fabs(realized_rewards_sum-best_reward_sum)
    realized_vs_best_ratio /= float(best_reward_sum)
    nose.tools.ok_( realized_vs_best_ratio < 0.25, 
                    'realized reward sum must be within 25% of best arm' )
