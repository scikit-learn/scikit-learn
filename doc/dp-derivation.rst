.. role:: raw-math(raw)
    :format: latex html

The spherical model
===================

 The model then is

$$
\begin{array}{rcl}
  \phi_k   &\sim& Beta(1, \alpha_1) \\
  \mu_k   &\sim& Normal(0,  \mathbf{I}) \\
  \sigma_k &\sim& Gamma(1, 1) \\
  z_{i}     &\sim& SBP(\phi) \\
  X_t &\sim& Normal(\mu_{z_i}, \frac{1}{\sigma_{z_i}} \mathbf{I})
\end{array}
$$

Tha variational distribution we'll use is
$$
\begin{array}{rcl}
  \phi_k   &\sim& Beta(\gamma_{k,1}, \gamma_{k,2}) \\
  \mu_k   &\sim& Normal(\nu_{\mu_k},  \mathbf{I}) \\
  \sigma_k &\sim& Gamma(a_{k}, b_{k}) \\
  z_{i}     &\sim& Discrete(\nu_{z_i}) \\
\end{array}
$$

The bound
--------

The variational bound is
$$
\begin{array}{rcl}
  \log P(X) &\ge& 
  \sum_k (E_q[\log P(\phi_k)] - E_q[\log Q(\phi_k)]) \\
  &&
  +\sum_k \left( E_q[\log P(\mu_k)] - E_q[\log Q(\mu_k)] \right) \\
  &&
  +\sum_k \left( E_q[\log P(\sigma_k)] - E_q[\log Q(\sigma_k)] \right) \\
  &&
  +\sum_i \left( E_q[\log P(z_i)] - E_q[\log Q(z_i)] \right) \\
  &&
  +\sum_i E_q[\log P(X_t)]
\end{array}
$$

**The bound for $\phi_k$**

$$
\begin{array}{rcl}
  E_q[\log Beta(1,\alpha)] - E[\log Beta(\gamma_{k,1},\gamma_{k,2})] 
  &=&
  \log \Gamma(1+\alpha) - \log \Gamma(\alpha) \\ && 
  +(\alpha-1)(\Psi(\gamma_{k,2})-\Psi(\gamma_{k,1}+\gamma_{k,2})) \\ &&
  - \log \Gamma(\gamma_{k,1}+\gamma_{k,2}) + \log \Gamma(\gamma_{k,1}) +
  \log \Gamma(\gamma_{k,2}) \\ &&
  -
  (\gamma_{k,1}-1)(\Psi(\gamma_{k,1})-\Psi(\gamma_{k,1}+\gamma_{k,2}))
  \\ &&
  -
  (\gamma_{k,2}-1)(\Psi(\gamma_{k,2})-\Psi(\gamma_{k,1}+\gamma_{k,2}))  
\end{array}
$$

**The bound for $\mu_k$**

$$
 \begin{array}{rcl}
  && E_q[\log P(\mu_k)] - E_q[\log Q(\mu_k)] \\
  &=&
  \int\!d\mu_f q(\mu_f) \log P(\mu_f)
    - \int\!d\mu_f q(\mu_f) \log Q(\mu_f)  \\
    &=&
- \frac{D}{2}\log 2\pi - \frac{1}{2} ||\nu_{\mu_k}||^2 - \frac{D}{2}
+ \frac{D}{2} \log 2\pi e 
\end{array}
$$


**The bound for $\sigma_k$**

Here I'll use the inverse scale parametrization of the gamma
distribution.

$$
\begin{array}{rcl}
  && E_q[\log P(\sigma_k)] - E_q [\log Q(\sigma_k)] \\ &=&
  \log \Gamma (a_k) - (a_k-1)\Psi(a_k) -\log b_k + a_k - \frac{a_k}{b_k}
\end{array}
$$


**The bound for z**

$$
\begin{array}{rcl}
  && E_q[\log P(z)] - E_q[\log Q(z)] \\
  &=&
  \sum_{k} \left( \left(\sum_{j=k+1}^K
      \nu_{z_{i,j}}\right)(\Psi(\gamma{k,1})-\Psi(\gamma{k,1}+\gamma_{k,2})) 
    +
    \nu_{z_{i,k}}(\Psi(\gamma_{k,1})-\Psi(\gamma_{k,1}+\gamma_{k,2}))\right)
  \\ &&
  - \sum_k
  \nu_{z_{i,k}} \log \nu_{z_{i,k}} \\
\end{array}
$$


**The bound for $X$**

Recall that there is no need for a $Q(X)$ so this bound is just

$$
\begin{array}{rcl}
  E_q[\log P(X_i)] &=& \sum_k \nu_{z_k} \left( - \frac{D}{2}\log 2\pi 
    +\frac{D}{2} (\Psi(a_k) - \log(b_k))
    -\frac{a_k}{2b_k} (||X_i - \nu_{\mu_k}||^2+D) - \log 2 \pi e  \right)
\end{array}
$$

For simplicity I'll later call the term inside the parenthesis $E_q[\log P(X_i|z_i=k)]$

The updates
----------

**Updating $\gamma$**

\begin{array}{rcl}
  \gamma_{k,1} &=& 1+\sum_i \nu_{z_{i,k}} \\
  \gamma_{k,2} &=& \alpha + \sum_i \sum_{j > k} \nu_{z_{i,j}}. 
\end{array}

**Updating $\mu$**

The updates for mu essentially are just weighted expectations of
$X$ regularized by the prior. We can see this by taking the
gradient of the bound w.r.t. $\nu_{\mu}$ and setting it to zero. The
gradient is

$$
  \nabla L = -\nu_{\mu_k} + \sum_i \frac{\nu_{z_{i,k}}b_k}{a_k}(X_i + -\nu_{\mu})
$$

so the update is
$$  \nu_{\mu_k} = \frac{\sum_i \frac{\nu_{z_{i,k}}b_k}{a_k}X_i}{1+\sum_i \frac{\nu_{z_{i,k}}b_k}{a_k}} $$


**Updating $a$ and $b$**


For some odd reason it doesn't really work when you derive the updates
for a and b using the gradients of the lower bound (it beats me why,
but some weird terms involving the $\Psi'$ function show up and it
gets really hard to isolate $a$). However, we can use the other formula,

$$\log Q(\sigma_k) = E_{v \ne \sigma_k}[\log P] + const $$

All the terms not involving $\sigma_k$ get folded over into the
constant and we get two terms: the prior and the probability of
$X$. This gives us

$$
\log Q(\sigma_k) = -\sigma_k  + \frac{D}{2} \sum_i \nu_{z_{i,k}}\log \sigma_k  - \frac{\sigma_k}{2}\sum_i \nu_{z_{i,k}} (||X_i-\mu_k||^2 + D)
$$

This is the log of a gamma distribution, with $a_k = 1+
\frac{D}{2}\sum_i \nu_{z_{i,k}}$ and

$$
b_k = 1 + \frac{1}{2}\sum_i \nu_{z_{i,k}} (||X_i-\mu_k||^2 + D).
$$


You can verify this by normalizing the previous term.

**Updating $z$**

$$\log \nu_{z_{i,k}} \propto \Psi(\gamma_{k,1}) -
\Psi(\gamma_{k,1} + \gamma_{k,2}) + E_Q[\log P(X_i|z_i=k)] +
\sum_{j < k} \left (\Psi(\gamma_{j,2}) -
  \Psi(\gamma_{j,1}+\gamma_{j,2})\right).
$$

The diagonal model
=================


The model then is

$$\begin{array}{rcl}
  \phi_k   &\sim& Beta(1, \alpha_1) \\
  \mu_k   &\sim& Normal(0,  \mathbf{I}) \\
  \sigma_{k,d} &\sim& Gamma(1, 1) \\
  z_{i}     &\sim& SBP(\phi) \\
  X_t &\sim& Normal(\mu_{z_i}, \bm{\sigma_{z_i}}^{-1})
\end{array}$$

Tha variational distribution we'll use is

$$\begin{array}{rcl}
  \phi_k   &\sim& Beta(\gamma_{k,1}, \gamma_{k,2}) \\
  \mu_k   &\sim& Normal(\nu_{\mu_k},  \mathbf{I}) \\
  \sigma_{k,d} &\sim& Gamma(a_{k,d}, b_{k,d}) \\
  z_{i}     &\sim& Discrete(\nu_{z_i}) \\
\end{array}
$$

The lower bound
--------------

The changes in this lower bound from the previous model are in the
distributions of $\sigma$ (as there are a lot more $\sigma$s now) and $X$.

The bound for $\sigma_{k,d}$ is the same bound for $\sigma_k$ and can
be safelly ommited.

**The bound for $X$**

The main difference here is that the precision matrix $\bm{\sigma_k}$
scales the norm, so we have an extra term after computing the
expectation of $\mu_k^T\bm{\sigma_k}\mu_k$, which is
$\nu_{\mu_k}^T\bm{\sigma_k}\nu_{\mu_k} + \sum_d \sigma_{k,d}$. We then
have

$$\begin{array}{rcl}
  E_q[\log P(X_i)] &=& \sum_k \nu_{z_k} \Big( - \frac{D}{2}\log 2\pi 
    +\frac{1}{2}\sum_d (\Psi(a_{k,d}) - \log(b_{k,d})) \\
    && 
    -\frac{1}{2}((X_i - \nu_{\mu_k})^T\bm{\frac{a_k}{b_k}}(X_i - \nu_{\mu_k})+ \sum_d \sigma_{k,d})- \log 2 \pi e  \Big)
\end{array}
$$

The updates
-----------

The updates only chance for $\mu$ (to weight them with the new
$\sigma$), $z$ (but the change is all folded into the
$E_q[P(X_i|z_i=k)]$ term), and the $a$ and $b$ variables themselves.

**The update for $\mu$**

$$  \nu_{\mu_k} = \left(\mathbf{I}+\sum_i \frac{\nu_{z_{i,k}}\mathbf{b_k}}{\mathbf{a_k}}\right)^{-1}\left(\sum_i \frac{\nu_{z_{i,k}}b_k}{a_k}X_i\right)
$$


**The updates for a and b**

Here we'll do something very similar to the spheric model. The main
difference is that now each $\sigma_{k,d}$ controls only one dimension
of the bound:

$$\log Q(\sigma_{k,d}) = -\sigma_{k,d} + \sum_i \nu_{z_{i,k}}\frac{1}{2}\log \sigma_{k,d} 
- \frac{\sigma_{k,d}}{2}\sum_i \nu_{z_{i,k}} ((X_{i,d}-\mu_{k,d})^2 + D)
$$

Hence 
$$ a_{k,d} = 1 + \frac{1}{2} \sum_i \nu_{z_{i,k}} $$
and
$$
b_{k,d} = 1 + \frac{1}{2} \sum_i \nu_{z_{i,k}}((X_{i,d}-\mu_{k,d})^2 + D).
$$

The tied model
=============

 The model then is
$$
\begin{array}{rcl}
  \phi_k   &\sim& Beta(1, \alpha_1) \\
  \mu_k   &\sim& Normal(0,  \mathbf{I}) \\
  \Sigma &\sim& Wishart(D, \mathbf{I}) \\
  z_{i}     &\sim& SBP(\phi) \\
  X_t &\sim& Normal(\mu_{z_i},  \Sigma^{-1})
\end{array}
$$
Tha variational distribution we'll use is

$$\begin{array}{rcl}
  \phi_k   &\sim& Beta(\gamma_{k,1}, \gamma_{k,2}) \\
  \mu_k   &\sim& Normal(\nu_{\mu_k},  \mathbf{I}) \\
  \Sigma &\sim& Wishart(a, \mathbf{B}) \\
  z_{i}     &\sim& Discrete(\nu_{z_i}) \\
\end{array}
$$

The lower bound
---------------

There are two changes in the lower-bound: for $\Sigma$ and for $X$.

**The bound for $\Sigma$**

$$
\begin{array}{rcl}
  \frac{D^2}{2}\log 2  + \sum_d \log \Gamma(\frac{D+1-d}{2}) \\
  - \frac{aD}{2}\log 2 + \frac{a}{2} \log |\mathbf{B}| + \sum_d \log \Gamma(\frac{a+1-d}{2}) \\
  + \frac{a-D}{2}\left(\sum_d \Psi\left(\frac{a+1-d}{2}\right) 
    + D \log 2 + \log |\mathbf{B}|\right) \\
  + \frac{1}{2} a \mathbf{tr}[\mathbf{B}-\mathbf{I}]
\end{array}
$$

**The bound for X**

$$
\begin{array}{rcl}
  E_q[\log P(X_i)] &=& \sum_k \nu_{z_k} \Big( - \frac{D}{2}\log 2\pi 
    +\frac{1}{2}\left(\sum_d \Psi\left(\frac{a+1-d}{2}\right) 
    + D \log 2 + \log |\mathbf{B}|\right) \\
    && 
    -\frac{1}{2}((X_i - \nu_{\mu_k})a\mathbf{B}(X_i - \nu_{\mu_k})+ a\mathbf{tr}(\mathbf{B}))- \log 2 \pi e  \Big)
\end{array}
$$

The updates
-----------

As in the last setting, what changes are the trivial update for $z$,
the update for $\mu$ and the update for $a$ and $\mathbf{B}$.

**The update for $\mu$**

$$  \nu_{\mu_k} = \left(\mathbf{I}+ a\mathbf{B}\sum_i \nu_{z_{i,k}}\right)^{-1}
    \left(a\mathbf{B}\sum_i \nu_{z_{i,k}} X_i\right)
$$

**The update for $a$ and $B$**

As this distribution is far too complicated I'm not even going to try
going at it the gradient way.

$$\log Q(\Sigma) = -\frac{1}{2}\log |\Sigma| - \frac{1}{2} \mathbf{tr}[\Sigma]
+ \sum_i \sum_k \nu_{z_{i,k}} \left( -\frac{1}{2}\log |\Sigma| - \frac{1}{2}((X_i-\nu_{\mu_k})^T\Sigma(X_i-\nu_{\mu_k})+\mathbf{tr}[I \Sigma]) \right)
$$

which non-trivially (seeing that the quadratic form with $\Sigma$ in
the middle can be expressed as the trace of something) reduces to

$$
\log Q(\Sigma) = -\frac{1}{2}\log |\Sigma| - \frac{1}{2} \mathbf{tr}[\Sigma]
+ \sum_i \sum_k \nu_{z_{i,k}} \left( -\frac{1}{2}\log |\Sigma| - \frac{1}{2}(\mathbf{tr}[(X_i-\nu_{\mu_k})(X_i-\nu_{\mu_k})^T\Sigma]+\mathbf{tr}[I \Sigma]) \right)
$$
hence this (with a bit of squinting) looks like a wishart with parameters

$$
a = 2 + D + T
$$
and
$$
\mathbf{B} = \left((T+1)\mathbf{I} + \sum_i \sum_k \nu_{z_{i,k}}(X_i-\nu_{\mu_k})(X_i-\nu_{\mu_k})^T\right)^{-1}
$$




The full model
=============

 The model then is

$$\begin{array}{rcl}
  \phi_k   &\sim& Beta(1, \alpha_1) \\
  \mu_k   &\sim& Normal(0,  \mathbf{I}) \\
  \Sigma_k &\sim& Wishart(D, \mathbf{I}) \\
  z_{i}     &\sim& SBP(\phi) \\
  X_t &\sim& Normal(\mu_{z_i},  \Sigma_{z,i}^{-1})
\end{array}
$$

Tha variational distribution we'll use is

$$\begin{array}{rcl}
  \phi_k   &\sim& Beta(\gamma_{k,1}, \gamma_{k,2}) \\
  \mu_k   &\sim& Normal(\nu_{\mu_k},  \mathbf{I}) \\
  \Sigma_k &\sim& Wishart(a_k, \mathbf{B_k}) \\
  z_{i}     &\sim& Discrete(\nu_{z_i}) \\
\end{array}
$$

The lower bound
--------------

All that changes in this lower bound in comparison to the previous one
is that there are K priors on different $\Sigma$ precision matrices
and there are the correct indices on the bound for X.

The updates
-----------

All that changes in the updates is that the update for mu uses only
the proper sigma and the updates for a and B don't have a sum over K, so 

$$  \nu_{\mu_k} = \left(\mathbf{I}+ a_k\mathbf{B_k}\sum_i \nu_{z_{i,k}}\right)^{-1}
    \left(a_k\mathbf{B_k}\sum_i \nu_{z_{i,k}} X_i\right)
$$

$$
a_k = 2 + D + \sum_i \nu_{z_{i,k}}
$$
and
$$
\mathbf{B} = \left(\left(\sum_i\nu_{z_{i,k}}+1\right)\mathbf{I} + \sum_i  \nu_{z_{i,k}}(X_i-\nu_{\mu_k})(X_i-\nu_{\mu_k})^T\right)^{-1}
$$

