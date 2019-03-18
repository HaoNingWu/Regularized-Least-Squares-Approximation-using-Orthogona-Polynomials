# Regularized-Least-Squares-Approximation-using-Orthogonal-Polynomials
This is a collection of MATLAB codes of numerical experiments in the paper 
"Regularized Weighted Discrete Least Squares Approximation by Orthogonal Polynomials" (by Congpei An and Haoning Wu)，
which is available on [arXiv:1805.01140](https://arxiv.org/abs/1805.01140).

To run these test codes, one should install [Chebfun](http://www.chebfun.org/), also can be ontained on [GitHub](https://github.com/chebfun), 
which is developed by Oxford's Numerical Analysis Group. See http://www.chebfun.org/ for Chebfun information.
We also provide Chenfun here! You are recommended to ADD IT TO PATH first.

Core codes:

1) l1_beta.m: solve l1-regularized approximation problem for coefficent beta
2) l2_beta.m: solve l2-regularized approximation problem for coefficent beta

Numerical experiments:

1) smoothapproximation.m：approximate a smooth function (This function stems from L. N. Trefethen's Approximation Theory and Approximation Practice)
2) reduction.m: recovery from a noisy function (spectral density of gate signal)
3) highoscillateddenoising.m: recovery from a highly oscillatory function with noise (function y=airy(40x))

Notes:

1) checkbary1.m: this code illustrates that it is the same numerically between l1 minimizer by closed-form solution and l1-regularized barycentric form
2) checkbary2.m: this code illustrates that it is the same numerically between l2 minimizer by closed-form solution and l2-regularized barycentric form
3) sparsity.m: this code shows the sparsity of l1-regularized solution

Contact information:
Hao-Ning Wu
Department of Mathematics, Jinan University 
Guangzhou, China 
haoning.wu@outlook.com, hwu@stu2015.jnu.edu.cn 
