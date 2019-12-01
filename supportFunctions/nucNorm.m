function [f,g]=nucNorm(X)

[u,d,v]=svd(X,'econ');
f=sum(diag(d));
g=u*v';
