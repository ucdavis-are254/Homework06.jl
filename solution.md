# Homework 06

In class, we discussed the L&L Example 6.4.1 growth model with and without an irreversible investment constraint: $f(k) \geq c$.

$$
\max_{c(t)} \int_0^T e^{-\rho t} u(c(t)) dt \qquad st \qquad \dot k = f(k) - c - mk
$$

Boundary conditions $k_0, k_T, T$

Assume $u(c) = \log c$ and $f(k) = k^\alpha$

Assume $\alpha = 0.5, \rho = 0.1, m = 0.3$.

Have $k_0 = 1.15$ and $k_T = 0.6$

Need to get a couple of $T$ values

# Solution

The current value Hamiltonian and Lagrangian are
$$
\tilde H(s,c,t) = u(c) + \psi[f(k) - c - mk]
$$

$$
\mathcal L = u(c) + \psi[f(k) - c - mk] + \lambda[f(k) - c]
$$

The FONC are
$$
\mathcal L_c = u'(c) - \psi - \lambda = 0
$$
The co-state is
$$
\dot \psi = [m+\rho]\psi -f'(k)[ \psi + \lambda]
$$
The state transition is
$$
\dot k = [f(k) - c] - mk
$$
It'll be convenient to write that
$$
u'(c) = \frac{1}{c} \qquad \qquad f'(k) = \alpha k^{-(1-\alpha)}
$$





## Unconstrained case

In the unconstrained case, we have $\lambda=0$, so
$$
u'(c) = \psi \qquad \dot \psi = [m+\rho - f'(k)]\psi \qquad \dot k = f(k) - c - mk
$$
Our convenient functional forms imply that $c = 1/\psi$, so we can write these 2 ODEs as
$$
\dot k = f(k) - 1/\psi - mk \qquad \qquad \dot \psi = [m+\rho - f'(k)]\psi
$$
The interior nullclines are
$$
\dot k \geq 0 \implies \psi \geq \frac{1}{f(k)-mk} \qquad \qquad \dot \psi \geq 0 \implies k \geq \left(\frac{\alpha}{m+\rho}\right)^\frac{1}{(1-\alpha)}
$$
The steady state is where these hold with equality


## Constrained case

The constraint boundary is where $c = f(k)$, which happens at $\psi = u'(f(k)) = \frac{1}{f(k)}$. At this point, we have that $\psi + \lambda = u'(f(k)) = \frac{1}{f(k)}$
$$
\dot k = -mk
\qquad \qquad
\dot \psi = (m + \rho)\psi - f'(k)u'(f(k))
$$
We have no constrained nullcline for $\dot k$ as this would require $k=0$. For $\psi$, we have that
$$
\dot \psi \geq 0 \implies \psi \geq \frac{f'(k)u'(f(k))}{(m+\rho)}
$$

### "scrap value" for constrained part

We also need to figure out the value of entering in the second part of the problem. Once we enter the constrained region, we have that
$$
\dot k = -mk \implies k(t) = k_T e^{m(T-t)}
$$
Suppose the time the constraint binds is $\hat t$. Then the integral is
$$
\int_{\hat t}^T e^{-\rho t} \log(k_T e^{m(T-t)}) dt
$$
We break this up into
$$
\int_{\hat t}^T e^{-\rho t} \log(k_T e^{m(T-t)}) dt = \log(k_Te^{mT}) \int_{\hat t}^T e^{-\rho t}dt - m \int_{\hat t}^T e^{-\rho t} t dt
$$
The first piece is
$$
\log(k_Te^{mT}) \int_{\hat t}^T e^{-\rho t}dt = \log(k_Te^{mT}) \frac{e^{-\rho \hat t} - e^{-\rho T}}{\rho}
$$
The second piece we use [integration by parts](https://en.wikipedia.org/wiki/Integration_by_parts)
$$
\int_a^b u(t)v'(t)dt = uv|_a^b - \int_a^b u'(t)v(t)dt
$$
<!-- |_ asdf -->
Take $u(t) = t$ and $v'(t) = e^{-\rho t} \implies v(t) = \frac{e^{-\rho t}}{-\rho}$. Then
$$
-m\int_{\hat t}^T e^{-\rho t} t dt
= -m\frac{e^{-\rho \hat t}(\rho \hat t + 1) - e^{-\rho T}(\rho T + 1)}{\rho ^2}
= -m \frac{e^{-\rho \hat t}\hat t - e^{-\rho T} T}{\rho} - \frac{m}{\rho}\frac{e^{-\rho \hat t}- e^{-\rho T}}{\rho}
$$
Thus, we have that the current value of being at $t=\hat t$ is this gross function:
$$
V(\hat t; T, k_T) = \int_{\hat t}^T e^{-\rho t} \log(k_T e^{m(T-t)}) dt  = \left[\log(k_T e^{mT}) - \frac m \rho \right]\frac{e^{-\rho \hat t}- e^{-\rho T}}{\rho} - m\frac{e^{-\rho \hat t}\hat t - e^{-\rho T} T}{\rho}
$$
We need to be able to compute $\frac{\partial V}{\partial \hat t}$ ~~and $\frac{\partial V}{\partial k(\hat t)}$~~ to get our TVC right. We can get
$$
\frac{\partial V}{\partial \hat t} = -\left[\log(k_T e^{mT}) - \frac m \rho \right]e^{-\rho \hat t} + \frac m \rho e^{-\rho \hat t}\left[\rho \hat t - 1\right] = e^{-\rho \hat t}\left[m \hat t - \log (k_Te^{mT}) - \frac{2m}{\rho} \right]
$$
~~We need $\partial V / \partial k_{\hat t}$. First, we know that $k(\hat t) = k_T e^{m(T-\hat t)}$, so $k_Te^{mT} = \hat k e^{m\hat t}$~~
$$
\partial V / \partial \hat k = \frac{\partial}{\partial \hat k} \left\{ \left[\log \hat k + m\hat t - \frac m \rho \right]\frac{e^{-\rho \hat t}- e^{-\rho T}}{\rho} - m\frac{e^{-\rho \hat t}\hat t - e^{-\rho T} T}{\rho} \right\} = \frac{1}{\hat k}\frac{e^{-\rho \hat t}- e^{-\rho T}}{\rho}
$$

### Free-time TVC for constrained part

The value for $\hat k = k(\hat t)$ is pinned down by the constraint: $k(\hat t) = \hat k = k_T e^{m(T-\hat t)}$

To pin down $\hat t$, we write this as a free-end time problem for $\hat t$. The present value TVC given in L&L 7.6 (p 245) is
$$
H(s,c,\pi,\hat t) + \frac{\partial \phi(b,\hat t)}{\partial \hat t} = 0.
$$
For us, we'll have $\phi(\hat t) = e^{-\rho \hat t} V(\hat t; T,k_T)$. So to convert this to current values, we multiply by $e^{\rho \hat t}$ and look at
$$
u(s,c) + \psi f(s,c) + \frac{ \partial V(\hat t;T,k_T)}{\partial \hat t} = 0
$$
We know that $\hat k \equiv k(\hat t)$. **Because of the continuity of $\psi$**, we can say that $\hat \psi = \psi(\hat t) = u'(f(\hat k))$. Now we have a cute root-finding problem for $\hat t$:
$$
\log f(\hat k) + \hat \psi[-m \hat k] + \frac{\partial V(\hat t; k_T, T)}{\partial \hat t} = 0
$$
where $\hat k = k_T e^{m(T-\hat t)}$ and $\hat \psi = u'(f(\hat k))$

### The BVP

We're now in a position to solve the BVP. We know $t_0,k_0$ and $\hat t, \hat k$. We don't know $\hat t$, but once we do, then we also know $\hat k, \hat \psi$. So, our problem is to find $\hat t$ such that $(\hat k - k_0) = \int_0^{\hat t} \dot k(t) dt$. Because we are searching over an unknown $\hat t$, this makes collocation tricky... so we'll try shooting.
