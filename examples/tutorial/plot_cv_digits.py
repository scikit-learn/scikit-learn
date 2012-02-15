import numpy as np
from scikits.learn import cross_val, datasets, svm

digits = datasets.load_digits()
X = digits.data
y = digits.target

svc = svm.SVC()
gammas = np.logspace(-6, -1, 10)

scores = list()
scores_std = list()
for gamma in gammas:
    svc.gamma = gamma
    this_scores = cross_val.cross_val_score(svc, X, y, n_jobs=-1)
    scores.append(np.mean(this_scores))
    scores_std.append(np.std(this_scores))

import pylab as pl
pl.figure(1, figsize=(2.5, 2))
pl.clf()
pl.axes([.1, .25, .8, .7])
pl.semilogx(gammas, scores)
pl.semilogx(gammas, np.array(scores) + np.array(scores_std), 'b--')
pl.semilogx(gammas, np.array(scores) - np.array(scores_std), 'b--')
pl.yticks(())
pl.ylabel('CV score')
pl.xlabel('gamma')
pl.ylim(0, 1.1)
#pl.axhline(np.max(scores), linestyle='--', color='.5')
pl.text(gammas[np.argmax(scores)], .9*np.max(scores), '%.3f' % np.max(scores), 
        verticalalignment='top',
        horizontalalignment='center',
        )


