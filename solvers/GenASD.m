function [ U,V , obj,time] = GenASD(M,opts )
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
    xTol=1e-6;
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
    fTol=1e-5;
end


bbopts.maxit = 50;
bbopts.xtol = 1e-5;
bbopts.gtol = 1e-3;
bbopts.ftol = 1e-5;
bbopts.alpha  = 1e-3;
bbopts.rho  = 1e-4;
bbopts.sigma  = 0.1;
bbopts.eta  = 0.8;


[m,n]=size(M);

U=rand(m,r);
V=rand(n,r);
U1=U;
V1=V;

P=[U;V];
[j] = find(M~=0);

b=M(j);

W=eye(r);
[row,col]=ind2sub([m,n], j);
tic;
for k=1:maxIter
    Pold=P;
    U2=U1;
    U1=U;
    
    V2=V1;
    V1=V;
    momentum=ex*min((k-3)/(k), 0);
    
    
    Ut=U+momentum*(U1-U2);
    
    UV=partXY(Ut', V', row, col, length(row));
    PX=sparse(row, col, UV'-b,m,n);
    d=2*Ut*W+mu*PX*V;      

    gV=partXY(d', V', row, col, length(row));
    
    t=(2*trace(d'*U*W)+(mu)*gV*(UV'-b))/(mu*norm(gV)^2+2*trace(d'*d*W));
    
    U=Ut-t*d;
    
    
    Vt=V+momentum*(V1-V2);

    UV=partXY(U', Vt', row, col, length(row));
    
    PX=sparse(row, col, UV'-b,m,n);

    d=2*Vt*W+ mu*PX'*U;

    Ug=partXY(U', d', row, col, length(row));
    
    
    t= (2*trace(d'*V*W)+(mu)*Ug*(UV'-b))/(mu*norm(Ug)^2+2*trace(d'*d*W));
    
    V=Vt-t*d;
    
    P=[U;V];
    
    [Vw,D]=eig(U'*U+V'*V);
    W=Vw*diag(f(diag(D)))*Vw';
    
    time(k)=toc;
    if isfield(opts, "obj")
        obj(k)=sum(diag(D)>0.01)+mu*norm(UV-b')^2;
        if k>1
            if abs(obj(k-1)-obj(k))/obj(k)<fTol
                break
            end
        end
    end
    if ((norm(P-Pold,'fro')/norm(Pold,'fro'))<xTol)
        break
    end
    
    
end
end


