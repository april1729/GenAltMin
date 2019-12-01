function [ P , obj] = gnrtn(P0,A,b,opts )
% Function Name: gncr_psd_exact
% This function solves the optimization problem:
%           min   sum_i f(sigma_i)
%           s.t.  A*vec(X)=b
%                 X psd
%
% Inputs:
%   P0 n by r matrix for initial variable
%   A constraint matrix
%   b constraint vector
%   opt options struct
%
% Outputs:
%   X optimal solution
%
% $Date: February 1, 2019
% ________________________________________

if isfield(opts, 'mu')
    mu=opts.mu;
else
    mu=10;
end
if isfield(opts, 'f')
    f=opts.f;
else
    f=@(x) 50^2./(50+x).^2;
end
if isfield(opts, 'maxIter')
    maxIter=opts.maxIter;
else
    maxIter=500;
end
if isfield(opts, 'xtol')
    xTol=opts.xtol;
else
    xTol=1e-10;
end
if isfield(opts, 'bbopts')
    bbopts=opts.bbopts;
else
    bbopts.maxit = 50;
    bbopts.xtol = 1e-4;
    bbopts.gtol = 1e-4;
    bbopts.ftol = 1e-4;
    bbopts.alpha  = 1e-3;
    bbopts.rho  = 1e-4;
    bbopts.sigma  = 0.1;
    bbopts.eta  = 0.8;
end
if isfield(opts, 'exact')
    exact=opts.exact;
else
    exact=true;
end


[n,r]=size(P0);

[~,nm]=size(A);


Lambda=zeros(size(b));
P=P0;
W=eye(r);
j = find(sum(A,1));
[rows,cols]=ind2sub([n,n], j);

for k=1:maxIter
    Pold=P;
    
    funIterk=@(x) objectivePSD(x,W,A,b,rows,cols,mu,n);
    %testGradient(P,funIterk)

    [P, obj2]= BBGradient(P, funIterk, bbopts);
    X=sparse_multiply(P,P, rows, cols, n,n);
    
    [V,D]=eig(P'*P);
    W=V*diag(f(diag(D)))*V';
    if exact
        Lambda=Lambda+(A*vec(X)-b);
    end
    
    
    if isfield(opts, 'obj')
        obj(k)=opts.obj(P);
    else
        obj(k)=norm(P-Pold,'fro')/norm(Pold,'fro');
    end
    if mod(k,10)==0
        disp("iteration "+k)
    end
    if ((norm(P-Pold,'fro')/norm(Pold,'fro'))<xTol)
        break
    end
end

X=P*P';

    function [f,g]=objectivePSD(P,W,A,b,rows,cols,mu,n)
        X=sparse_multiply(P,P, rows, cols, n,n);
        f=trace(W*(P'*P))+(mu/2)*norm(A*vec(X)-b)^2;
        g=2*P*W+ 2*mu*reshape(A'*(A*vec(X)-b),[n,n])*P;
    end

end

