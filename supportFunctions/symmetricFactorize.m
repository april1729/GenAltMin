function [V]=symmetricFactorize(X,r)

[U,D]=eigs(X,r);

V=U*sqrt(D);