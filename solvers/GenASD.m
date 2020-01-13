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
    xTol=0;
end
if isfield(opts, 'gamma')
    gamma=opts.gamma;
else
    gamma=1000;
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
    fTol=1e-8;
end

if isfield(opts, 'scaled')
    scaled=opts.scaled;
else
    scaled=0;
end

[m,n]=size(M);

U=rand(m,r);
V=rand(n,r);
U1=U;
V1=V;

P=[U;V];
[j] = find(M~=0);

b=M(j);
p=length(b)/(m*n);
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
    if scaled 
        d=d*inv((1/p)*mu*V'*V+W);
    end
    gV=partXY(d', V', row, col, length(row));
    
    t=(2*trace(d'*U*W)+(mu)*gV*(UV'-b))/(mu*norm(gV)^2+2*trace(d'*d*W));
    
    U=Ut-t*d;
    
    
    Vt=V+momentum*(V1-V2);

    UV=partXY(U', Vt', row, col, length(row));
    
    PX=sparse(row, col, UV'-b,m,n);

    d=2*Vt*W+ mu*PX'*U;
    if scaled 
        d=d*inv((1/p)*mu*U'*U+W);
    end

    Ug=partXY(U', d', row, col, length(row));
    
    
    t= (2*trace(d'*V*W)+(mu)*Ug*(UV'-b))/(mu*norm(Ug)^2+2*trace(d'*d*W));
    
    V=Vt-t*d;
    
    P=[U;V];
    gamma=max(gamma*0.95,0.1);
    [Vw,D]=eig(U'*U+V'*V);
    W=Vw*diag(f(diag(D), gamma))*Vw';
    time(k)=toc;
    if isfield(opts, "obj")
        obj(k)=opts.obj(U,V);
        if k>1
            abs(obj(k-1)-obj(k))/obj(k);
            if (obj(k-1)-obj(k))/obj(k)<fTol
                break
            end
        end
    end
    if ((norm(P-Pold,'fro')/norm(Pold,'fro'))<xTol)
        break
    end
    
    
end
end


