function [ X , obj] = gnrtn(P0,A,b,opts )
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
    xTol=1e-8;
end
if isfield(opts, 'bbopts')
    bbopts=opts.bbopts;
else
    bbopts.maxit = 1000;
    bbopts.xtol = 1e-8;
    bbopts.gtol = 1e-8;
    bbopts.ftol = 1e-10;
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
if nm==n^2 && n>r
    isPSD=true;
else
    isPSD=false;
end


Lambda=zeros(size(b));
P=P0;
W=eye(r);
[~,j] = find(A);
[rows,cols]=ind2sub([n,n], j);

for k=1:maxIter
    Pold=P;
    if isPSD
        funIterk=@(x) objectivePSD(x,W,A,b,rows,cols,mu,n);
    else
        funIterk=@(x) objectiveRect(x,W,A,b,mu,n,r);
    end
    
    [P, obj2]= BBGradient(P, funIterk, bbopts);
    if isPSD
        X=P*P';
    else
        X=P;
    end
    [V,D]=eig(P'*P);
    W=V*diag(f(diag(D)))*V';
    if exact
        Lambda=Lambda+(A*vec(X)-b);
    end
    obj(k)=sum(diag(D)>0.001);
    
    if ((norm(P-Pold,'fro')/norm(Pold,'fro'))<xTol)
        break
    end
end

    function [f,g]=objectivePSD(P,W,A,b,rows,cols,mu,n)
        f=trace(W*(P'*P))+(mu/2)*norm(A*vec(P*P')-b,'fro')^2;
        g=2*P*W+ 2*mu*reshape(A'*(A*vec(P*P')-b),[n,n])*P;
    end

    function [f,g]=objectiveRect(X,W,A,b,mu,n,m)
        try
            f=trace(W*(X'*X))+(mu/2)*norm(A*vec(X)-b,'fro')^2;
            g=2*X*W+ 2*mu*reshape(A'*(A*vec(X)-b),[n,m]);
        catch
            disp('oh no')
        end
    end

end

