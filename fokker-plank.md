Writing a periodic progress report is a good way to record the process that you gain knowledge. Reviewing my progress reports in the last 3 years, it recovers my memory about the knowledge and technologies that I almost forgot as I rarely use them. From this post on, I would like to give them a review and share on the blog. Hope it can help people who are interest in and one can help me if I am wrong at somewhere.

Derivation
Here I briefly present the derivation for the Fokker-Planck equation from a stochastic differential equation.

Given the stochastic process

$$
dx=a(x,t)dt+b(x,t)dW_{t}
$$

where \(W_{t}\) is a Wiener process. By Ito lemma, for any twice-differentiable scalar function \(f(x)\) we have

$$
df(x)=\left(a(x,t)f’(x)+\frac{1}{2}b^{2}(x,t)f’’(x)\right)dt+b(x,t)f’(x)dW_{t}
$$

The expectation of \(f(x,t)\) yields

$$
E(f(x))=\int f(x)p(x,t)dx
$$

and take the derivative

$$
\frac{dE(f(x))}{dt}=\frac{d\int f(x)p(x,t)dx}{dt}=\int f(x)\frac{\partial p(x,t)}{\partial t}dx \tag{1}
$$

Also, we could plug Eq.[1] in the expectation of \(f(x)\) and take the derivative yields

$$
\frac{dE(f(x))}{dt}=\frac{E(df(x))}{dt}=E\left(a(x,t)f’(x)+\frac{1}{2}b^{2}(x,t)f’’(x)\right) \tag{2}
$$

From that Eq.[1] and Eq.[2] are identical, we have

$$
\begin{align}
\int f(x)\frac{\partial p(x,t)}{\partial t}dx &=\int\left(a(x,t)f’(x)+\frac{1}{2}b^{2}(x,t)f’’(x)\right)p(x,t)dx \\
& =\int a(x,t)f’(x)p(x,t)dx+\frac{1}{2}\int b^{2}(x,t)f’’(x)p(x,t)dx \\
& =-\int f(x)\frac{\partial a(x,t)p(x,t)}{\partial x}dx+\frac{1}{2}\int f(x)\frac{\partial^{2}b^{2}(x,t)p(x,t)}{\partial x^{2}}dx \\
& =\int f(x)\left(-\frac{\partial a(x,t)p(x,t)}{\partial x}+\frac{1}{2}\frac{\partial^{2}b^{2}(x,t)p(x,t)}{\partial x^{2}}\right)dx
\end{align}
$$

As \(f(x)\) is arbitrary, we obtain the Fokker-Planck equation in one dimension

$$
\frac{\partial p(x,t)}{\partial t}=-\frac{\partial a(x,t)p(x,t)}{\partial x}+\frac{1}{2}\frac{\partial^{2}b^{2}(x,t)p(x,t)}{\partial x^{2}}
$$
