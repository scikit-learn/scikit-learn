Insights for Enet dual gap computation in `sklearn/linear_model`
================================================================

Let us start by reminding duality properties (mostly extracted from [1]) of the
Lasso problems, with `sklearn` notation:

- :math:`n = n_samples`
- :math:`p = n_features`
- :math:`X = [x_1, \dots, x_p]`

Lasso reminder
--------------

.. math::
   :nowrap:

   \mbox{Primal Lasso problem: } \quad
   w^{(\alpha,X,y)} \in \argmin_{w \in \bbR^p}\underbrace{\frac{1}{2} \norm{X
   w - y}_2^2 + \alpha \normin{w}_1}_{=\cP^{(\alpha,X,y)}(w)} \enspace .


Denoting :math:`\Delta_X = \big\{ \theta \in \bbR^n \; : \;
\norm{X^\top \theta }_\infty \leq 1\big\}\subset \bbR^n` the dual feasible set,
a dual formulation of the Lasso reads (see for instance [2] or [3]):

.. math::
   :nowrap:

   \mbox{Dual Lasso problem: }
   \theta^{(\alpha,X,y)}=\argmax_{\theta \in \Delta_X}
   \underbrace{\frac{1}{2}\norm{y}^2_2 - \frac{\alpha^2}{2}\norm{\theta -
   \frac{y}{\alpha}}^2_2}_{=\cD^{(\alpha,y)}(\theta)}.

With the residual notation :math`r=y-X w`, and :math:`\theta = r/c` (:math`c`
is a constant determined later) the duality gap reads:

.. math::
   :nowrap:

   \begin{align}
   \mbox{Duality gap: } \quad
   \cG^{(\alpha,X,y)}(w,\theta)=&\cP^{(\alpha,X,y)}(w)-\cD^{(\alpha,y)}(\theta)\\
   =&\frac{1}{2} \norm{X w - y}_2^2 + \alpha \norm{w}_1
   -\frac{1}{2}\norm{y}^2_2 + \frac{\alpha^2}{2}\norm{\theta -
   \frac{y}{\alpha}}^2_2\\
   =&\frac{1}{2} \norm{r}_2^2 + \alpha \norm{w}_1
   + \frac{\alpha^2}{2}\norm{\theta}^2_2 - \alpha\theta^\top y\\
   =&\frac{1}{2} \norm{r}_2^2 \left(1 + \frac{\alpha^2}{c^2}\right) + \alpha \norm{w}_1- \frac{\alpha}{c} r^\top y
   \enspace.
   \end{align}

Dual point creation: we choose a dual point :math:`\theta` that is proportional
to the residual :math:`r` (since at optimality :math:`\theta^{(\alpha,X,y)}=(y-Xw^{(\alpha,X,y)})/\alpha`)

.. math::
   :nowrap:

   \theta=\frac{r}{c}, \text{ where } c=\max\left(\norm{X^\top r}_\infty,\alpha\right) \enspace.

this means we normalize :math`r/\alpha` only if this vector is not dual feasible
(this avoid rescaling issues when :math`X^\top r` is to close from zero).

Positive Lasso
--------------
In the case where positive constraints are added to the Lasso,
the primal problem reads

.. math::
   :nowrap:

   \mbox{Primal Positive Lasso problem: } \quad
   w^{(\alpha,X,y)}_P \in \argmin_{\substack{w \in \bbR^p \\  \forall j \in [p], w_j \geq 0}}
   \frac{1}{2} \norm{X w - y}_2^2 + \alpha \norm{w}_1 \enspace .

.. math::
   :nowrap:

   \mbox{Dual Lasso problem: }
   \theta^{(\alpha,X,y)}_P, \mu^{(\alpha,X,y)}_P=\argmax_{\substack{\theta \in \bbR^n, \mu \in \bbR^p\\
   \norm{X^\top \theta - \mu }_\infty \leq 1}\\
   \forall j \in [p], \mu_j\leq 0}
   \frac{1}{2}\norm{y}^2_2 - \frac{\alpha^2}{2}\norm{\theta - \frac{y}{\alpha}}^2_2.

Hence, the same dual point can be considered for $\theta$ as for the Lasso,
and the simple choice :math`\mu=0`, leads to a dual feasible couple, and duality
gap evaluations then coincide.

From Lasso to enet
------------------

Let us write enet problem:

.. math::
   :nowrap:

   \mbox{Primal Enet problem: } \quad
   w^{(\alpha,\beta,X,y)}_{enet} \in \argmin_{w \in \bbR^p}
   \underbrace{\frac{1}{2} \norm{X w - y}_2^2 + \alpha \norm{w}_1 +\beta \frac{\norm{w}^2}{2}}_{\cP_{enet}^{(\alpha,\beta,X,y)}(w)}
   \enspace .

Let us rewrite the enet problem as a Lasso formulation.

For that we need:
:math`\tilde{X}=
\begin{pmatrix}
X       \\
\sqrt{\beta} \Id_{p}
\end{pmatrix} \in \bbR^{(n+p) \times p}`
,
:math`\tilde{y}=
\begin{pmatrix}
y       \\
0_{p}
\end{pmatrix} \in \bbR^{n+p}`,
and
:math`\tilde{r}=
\begin{pmatrix}
r       \\
-\sqrt{\beta} w
\end{pmatrix} \in \bbR^{n+p}`,

So now one can notice that

.. math::
   :nowrap:

   \cP_{enet}^{(\alpha,\beta,X,y)}(w)=\cP^{(\alpha,\tilde{X},\tilde{y})}(w)

At optimality we can defined the dual optimal point:

.. math::
   :nowrap:

   \begin{align}
   \tilde{\theta}^{(\alpha,\beta,X,y)}=&\frac{\tilde{y}-X\tilde{w}^{(\alpha,\beta,X,y)}}{\alpha}\\
   =& \begin{pmatrix}
   \frac{\tilde{y}-X w^{(\alpha,X,y)}}{\alpha}\\
   -\frac{\sqrt{\beta}}{\alpha}w^{(\alpha,X,y)}
   \end{pmatrix}
   \end{align}

Hence we propose as a dual feasible point
:math`\tilde{\theta}=\frac{\tilde{r}}{\tilde{c}}`` where we define


.. math::
   :nowrap:
   \tilde{c}=\max\left(\norm{\tilde{X}^\top \tilde{r}}_\infty,\alpha\right)=\max\left(\norm{X^\top r-\beta w}_\infty,\alpha\right)\enspace.

meaning

.. math::
   :nowrap:

   \tilde{\theta}=
   \frac{1}{\tilde{c}}\begin{pmatrix}
   r\\
   -\sqrt{\beta}w
   \end{pmatrix}

Hence we can write the enet duality gap as:

.. math::
   :nowrap:
   \begin{align}
   \cG_{enet}^{(\alpha,\beta,X,y)}(\beta, \theta)&=
   \cP^{(\alpha,\tilde{X},\tilde{y})}(w)-\cD^{(\alpha,\tilde{y})}(\tilde{\theta})
   \\
   &=\frac{1}{2} \norm{\tilde{r}}_2^2 \left(1 + \frac{\alpha^2}{\tilde{c}^2}\right) + \alpha \norm{w}_1- \frac{\alpha}{\tilde{c}} \tilde{r}^\top {\tilde{y}}\\
   \cG_{enet}^{(\alpha,\beta,X,y)}(\beta, \theta)
   &=\frac{1}{2} \left(\norm{r}_2^2+\beta \norm{w}^2\right) \left(1 + \frac{\alpha^2}{\tilde{c}^2}\right) + \alpha \norm{w}_1- \frac{\alpha}{\tilde{c}} r^\top y
   \enspace.
   \end{align}
