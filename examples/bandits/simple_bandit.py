"""
Simple example of working bandit application
"""

from sklearn.bandits.context_free_bandits import EpsilonGreedyBandit

# let's generate reward data for 3 possible actions
arms_rewards = [ 
                [0 for _ in range(0,10)], # arm 0 performs worst
                [1 for _ in range(0,10)], # arm 1 performs best
                [0,1,0,1,0,1,0,1,0,1]     # arm 2 is mixed
                ]

bandit = EpsilonGreedyBandit()

# let's assume a single user
visitor='lone'
arms = [0,1,2]
realized_reward = 0

for _ in range(0,10):
    chosen_action = bandit.recommend(visitor, 
                                     arms)
    reward = arms_rewards[chosen_action][_]
    realized_reward += reward
    bandit.update(visitor, chosen_action, reward)

print "realized reward sum:", realized_reward, "/ max. possible:", 10 