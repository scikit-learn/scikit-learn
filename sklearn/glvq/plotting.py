import matplotlib.pyplot as plt

def project_plot2d(model, X, y,figure,title=""):
    dim = 2
    f = plt.figure(figure)
    f.suptitle(title)
    pred = model.predict(X)

    if hasattr(model,'psis_'):
        nb_prototype = model.w_.shape[0]
        ax = f.add_subplot(1,nb_prototype+1,1)
        ax.scatter(X[:, 0], X[:, 1], c=y, alpha=0.5)
        ax.scatter(X[:, 0], X[:, 1], c=pred, marker='.')
        ax.scatter(model.w_[:, 0], model.w_[:, 1])
        ax.axis('equal')
        for i in range(nb_prototype):
            X_p = model.project(X, i, dim)
            w_p = model.project(model.w_[i], i, dim)

            ax = f.add_subplot(1,nb_prototype+1,i+2)
            ax.scatter(X_p[:, 0], X_p[:, 1], c=y, alpha=0.2)
            #ax.scatter(X_p[:, 0], X_p[:, 1], c=pred, marker='.')
            ax.scatter(w_p[0], w_p[1],marker='D',s=20)
            ax.axis('equal')

    else:
        ax = f.add_subplot(121)
        ax.scatter(X[:, 0], X[:, 1], c=y, alpha=0.5)
        ax.scatter(X[:, 0], X[:, 1], c=pred, marker='.')
        ax.scatter(model.w_[:, 0], model.w_[:, 1])
        ax.axis('equal')
        X_p = model.project(X, dim)
        w_p = model.project(model.w_, dim)

        ax = f.add_subplot(122)
        ax.scatter(X_p[:, 0], X_p[:, 1], c=y, alpha=0.5)
        #ax.scatter(X_p[:, 0], X_p[:, 1], c=pred, marker='.')
        ax.scatter(w_p[:, 0], w_p[:, 1],marker='D',s=20)
        ax.axis('equal')
    f.show()
