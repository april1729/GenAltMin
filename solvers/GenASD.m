function [ U,V , obj,time] = GenASD(X0,A,b,opts )
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
%_____________________________________________

obj=[];time=[];
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
    maxIter=10000;
end
if isfield(opts, 'xTol')
    xTol=opts.xTol;
else
    xTol=1e-5;
end
if isfield(opts, 'SDopts')
    SDopts=opts.SDopts;
else
    SDopts.maxIter=100;
    SDopts.xTol=1e-5;
    SDopts.fTol=1e-5;
end

if isfield(opts, 'exact')
    exact=opts.exact;
else
    exact=false;
end

if isfield(opts, 'r')
    r=opts.r;
else
    r=10;
end

if isfield(opts, 'ex')
    ex=opts.ex;
else
    ex=0;
end

if isfield(opts, 'fTol')
    fTol=opts.fTol;
else
    fTol=1e-4;
end


bbopts.maxit = 50;
bbopts.xtol = 1e-5;
bbopts.gtol = 1e-3;
bbopts.ftol = 1e-5;
bbopts.alpha  = 1e-3;
bbopts.rho  = 1e-4;
bbopts.sigma  = 0.1;
bbopts.eta  = 0.8;


[m,n]=size(X0);

[U, D, V]=svds(X0, r);
U=U*sqrt(D);
V=V*sqrt(D);
U1=U;
V1=V;

P=[U;V];

Lambda=zeros(size(b));
W=eye(r);
[~,j] = find(A);
[rows,cols]=ind2sub([m,n], j);
tic;
for k=1:maxIter
    Pold=P;
    U2=U1;
    U1=U;
    
    V2=V1;
    V1=V;
    momentum=ex*min((k-3)/(k), 0);
    
    
    Ut=U+momentum*(U1-U2);
    [~,d]= objectiveU(Ut,V,W,A,b,mu);
    
    %testGradient(U,@(x) objectiveU(x,V,W,A,b,mu));
    t=tu(Ut,V,W,A,b,mu,d);
    
    U=Ut-t*d;
    
    
    Vt=V+momentum*(V1-V2);
    [~,d] = objectiveV(U,Vt,W,A,b,mu);
    %testGradient(V,@(x) objectiveV(U,x,W,A,b,mu));
    
    t= tv(U,Vt,W,A,b,mu,d);
    
    V=Vt-t*d;
    
    P=[U;V];
    
    [Vw,D]=eig(U'*U+V'*V);
    W=Vw*diag(f(diag(D)))*Vw';
    X=sparse_multiply(U,V,rows, cols, m,n);
    
    time(k)=toc;
    if isfield(opts, "obj")
        obj(k)=opts.obj(U,V);
        if k>1
            if (obj(k-1)-obj(k))/obj(k)<fTol
                break
            end
        end
    end
    if ((norm(P-Pold,'fro')/norm(Pold,'fro'))<xTol)
        break
    end
    
    
end


    function [f,g]=objectiveU(x,V,W,A,b,mu)
        f=trace(W*(x'*x))+(mu/2)*norm(A*vec(sparse_multiply(x,V, rows, cols,m,n))-b,'fro')^2;
        g=2*x*W+mu*reshape(A'*(A*vec(sparse_multiply(x,V, rows, cols,m,n))-b),[m,n])*V;
    end

    function [f,g]=objectiveV(U,x,W,A,b,mu)
        f=trace(W*(x'*x))+(mu/2)*norm(A*vec(sparse_multiply(U,x, rows, cols,m,n))-b,'fro')^2;
        g=2*x*W+ mu*reshape(A'*(A*vec(sparse_multiply(U,x, rows, cols,m,n))-b),[m,n])'*U;
    end


    function [t]=tu(U,V,W,A,b,mu,g)
        UV=sparse_multiply(U,V, rows, cols,m,n);
        gV=sparse_multiply(g,V, rows, cols,m,n);
        t=(2*trace(g'*U*W)+(mu)*(A*vec(gV))'*(A*(vec(UV))-b))/(mu*norm(A*vec(gV), 'fro')^2+2*trace(g'*g*W));
    end

    function t=tv(U,V,W,A,b,mu,g)
        
        
        UV=sparse_multiply(U,V, rows, cols,m,n);
        Ug=sparse_multiply(U,g, rows, cols,m,n);
        t=(2*trace(g'*V*W)+(mu)*(A*vec(Ug))'*(A*(vec(UV))-b))/(mu*norm(A*vec(Ug), 'fro')^2+2*trace(g'*g*W));
        
    end
end


