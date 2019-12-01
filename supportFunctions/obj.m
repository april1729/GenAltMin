function [f_out,g_out]=obj(x, f,g, mu, A, b)

[u,d,v]=svd(x, 'econ');

g_out=u*diag(g(diag(d)))*v'+mu*reshape(A'*(A*vec(x)-b),size(x));
f_out=sum(f(diag(d)))+(mu/2)*norm(A*vec(x)-b)^2;