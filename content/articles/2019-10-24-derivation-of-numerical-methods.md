Title: Derivation of Numerical Methods
date: 2019-10-24 10:00
Modified: 2019-10-24 22:12
comments: true
slug: cpp-101
tags: implicit euler, explicit euler

<!-- PELICAN_BEGIN_SUMMARY -->
In this section, we summarize the second session about Derivation of Numerical Method. Let's go!
<!-- PELICAN_END_SUMMARY -->

## Some Elementary Schemes

* We consider first-order systems of ODEs in explicit form

$$
u^{\prime}=f\left(t,u\left(t\right)\right),\quad
t\in\left(t_{0},t_{0}+T\right],\quad
u\left(t_{0}\right)=u_{0}\tag{11}
$$

to determine the unknown function $u\colon\left[t_{0},t_{0}+T\right]\rightarrow\mathbb{R}^{d}$.

* In components this reads:

$$
\begin{pmatrix}
u^{\prime}_{1}\left(t\right)\\
\vdots\\
u^{\prime}_{d}\left(t\right)
\end{pmatrix}=
\begin{pmatrix}
f_{1}\left(t,u_{1}\left(t\right),\ldots,u_{d}\left(t\right)\right)\\
\vdots\\
f_{d}\left(t,u_{1}\left(t\right),\ldots,u_{d}\left(t\right)\right)
\end{pmatrix}.
$$

* The right hand side $f\colon\left[t_{0},t_{0}+T\right]\times\mathbb{R}^{d}\rightarrow\mathbb{R}^d$ is Lipschitz-continuous

$$
\left\|f\left(t,u\right)-f\left(t,w\right)\right\|\leq L\left(t\right)\left\|u-w\right\|.
$$

* Thus, the system has a unique solution.

## Taylor's Theorem

An important tool in the derivation and analysis of numerical methods for ODEs is the following theorem.

### Theorem (Taylor's Theorem with Lagrangian Remainder)

Let $u\colon I\rightarrow\mathbb{R}$ be $\left(n+1\right)$-times continuosly differentiable. Then, for $t$, $t+\Delta t\in I$, it holds

$$
u\left(t+\Delta t\right)=
\sum_{k=0}^{n}\frac{u^{\left(k\right)}\left(t\right)}{k!}\Delta t^{k}+\frac{u^{\left(n+1\right)}\left(t+\theta\Delta t\right)}{\left(n+1\right)!}\Delta t^{n+1},\quad\theta\in\left[0,1\right].
$$

As a consequence of Taylor's theorem we have for $n=1$.

$$
u\left(t+\Delta t\right)
=u\left(t\right)+u^{\prime}\left(t\right)\Delta t+\frac{u^{\prime\prime}\left(t+\xi\right)}{2}\Delta t^{2},\quad
0\leq\xi\leq\Delta t.
$$

Taken component-wise this hold for vector-valued $u$.

### Explicit Euler Method

* Choose $N$ time steps

$$
t_{0}<t_{1}<t_{2}<\cdots<t_{N-1}<t_{N}
=t_{0}+T,\quad\Delta t_{i}=t_{i+1}-t_{i}.
$$

* $y^{N}_{i}$ denotes the approximation of $u\left(t_{i}\right)$ computed with $N$ steps.

* Take Taylor, use ODE and omit error term to obtain the explicit Euler approximation

$$
\begin{aligned}
u\left(t_{i+1}\right)
&=u\left(t_{i}\right)+\Delta t_{i}u^{\prime}\left(t_{i}\right)+\frac{t_{i}+\xi_{i}}{2}\Delta t^{2}_{i}\\
&=u\left(t_{i}\right)+\Delta t_{i}f\left(t_{i},u\left(t_{i}\right)\right)+\frac{t_{i}+\xi_{i}}{2}\Delta t^{2}_{i}\\
\implies y^{N}_{i+1}
&=y^{N}_{i}+\Delta t_{i}f\left(t_{i},y^{N}_{i}\right).
\end{aligned}
$$

* Assuming $y^{N}_{i}=u\left(t_{i}\right)$ and substracting we obtain

$$
u\left(t_{i+1}\right)-y^{N}_{i+1}=\frac{u^{\prime\prime}\left(t_{i}+\xi_{i}\right)}{2}\Delta t^{2}_{i}
$$

* The error after one step is $O\left(\Delta t^2\right)$, how does it propagate?

### Implicit Euler Method

* Using Taylor's theorem slightly differently gives

$$
\begin{aligned}
u\left(t_{i}\right)
&=u\left(t_{i+1}-\Delta t_{i}\right)=u\left(t_{i+1}\right)-\Delta t_{i}u^{\prime}\left(t_{i+1}\right)+\Delta t^{2}_{i}\frac{u^{\prime\prime}\left(t_{i+1}-\xi_{i}\right)}{2}\\
&=u\left(t_{i+1}\right)-\Delta t_{i}f\left(t_{i+1},u\left(t_{i+1}\right)\right)+\Delta t^{2}_{i}\frac{u^{\prime\prime}\left(t_{i+1}-\xi_{i}\right)}{2}\\
\iff u\left(t_{i+1}\right)&-\Delta t_{i}f\left(t_{i+1},u\left(t_{i+1}\right)\right)=u\left(t_{i}\right)-\Delta t^{2}_{i}\frac{u^{\prime\prime}\left(t_{i+1}-\xi_{i}\right)}{2}
\end{aligned}
$$

* Which yields the implicit Euler approximation

$$
y^{N}_{i+1}-\Delta t_{i}f\left(t_{i+1},y^{N}_{i+1}\right)=y^{N}_{i}\tag{12}
$$

* Need to solve a nonlinear algebraic equation to obtain $y^{N}_{i+1}$ which is computationally much more demanding!


### Local Error in Implicit Euler Method

* We can modify the analysis of the explicit scheme.

* From the construction of the scheme we obtain.

$$
\begin{aligned}
u\left(t_{i+1}\right)
&=u\left(t_{i}\right)+\Delta t_{i}f\left(t_{i+1}, u\left(t_{i+1}\right)\right)-\Delta t^{2}_{i}\frac{u^{\prime\prime}\left(t_{i+1}-\xi_{i}\right)}{2}\\
y^{N}_{i+1}
&=y^{N}_{i}+\Delta t_{i}f\left(t_{i},y^{N}_{i+1}\right)
\end{aligned}
$$

* Substracting, taking norms and using $L$-continuity gives

$$
\begin{aligned}
u\left(t_{i+1}\right)-y^{N}_{i+1}
&=u\left(t_{i}\right)-y^{N}_{i}+\Delta t_{i}\left[f\left(t_{i+1},u\left(t_{i+1}\right)\right)-f\left(t_{i},y^{N}_{i+1}\right)\right]-\Delta t^{2}_{i}\frac{u^{\prime\prime}\left(t_{i+1}-\xi_{i}\right)}{2}\\
\left\|u\left(t_{i+1}\right)-y^{N}_{i+1}\right\|
&\leq\left\|u\left(t_{i}\right)-y^{N}_{i}\right\|+\Delta t_{i}L\left(t_{i+1}\right)\left\|u\left(t_{i+1}\right)-y^{N}_{i+1}\right\|+\frac{\Delta t^{2}_{i}}{2}\left\|u^{\prime\prime}\left(t_{i+1}-\xi_{i}\right)\right\|
\end{aligned}
$$

* For $\Delta t<\frac{1}{L}$ the error after one step is also $O\left(\Delta t^2\right)$.

* This time step restriction is actually superficial.

### Implicit Trapezoidal Rule

* Another approach follows from integrating the ODE with the trapezoidal rule

$$
u\left(t_{i+1}\right)-u\left(t_{i}\right)=\int\limits_{t_{i}}^{t_{i+1}}f\left(\xi,u\left(\xi\right)\right)d\xi=\frac{\Delta t_{i}}{2}\left[f\left(t_{i},u\left(t_{i}\right)\right)+f\left(t_{i+1},u\left(t_{i+1}\right)\right)\right]+O\left(\Delta t^{3}_{i}\right).
$$

* Resulting in the *implicit trapezoidal rule*

$$
\begin{aligned}
y^{N}_{i+1}-y^{N}_{i}
&=\frac{\Delta t_{i}}{2}\left[f\left(t_{i},y^{N}_{i}\right)+f\left(t_{i+1},y^{N}_{i+1}\right)\right]\\
\iff y^{N}_{i+1}-\frac{\Delta t_{i}}{2}.
\end{aligned}
$$