function [ U,V , obj,time] = GenAltMin(X0,A,b,opts )
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
if isfield(opts, 'beta')
    beta=opts.beta;
else
    beta=10;
end
if isfield(opts, 'f')
    f=opts.f;
else
    f=@(x,gamma) gamma^2./(gamma+x).^2;
end
if isfield(opts, 'maxIter')
    maxIter=opts.maxIter;
else
    maxIter=100;
end
if isfield(opts, 'xTol')
    xTol=opts.xTol;
else
    xTol=1e-5;
end
if isfield(opts, 'SDopts')
    SDopts=opts.SDopts;
else
    SDopts.maxIter=25;
    SDopts.xTol=1e-5;
    SDopts.fTol=1e-5;
end


if isfield(opts, 'r')
    r=opts.r;
else
    r=10;
end


if isfield(opts, 'fTol')
    fTol=opts.fTol;
else
    fTol=1e-5;
end
if isfield(opts, 'gamma')
    gamma=opts.gamma;
else
    gamma=100;
end

[m,n]=size(X0);

[U, D, V]=svds(X0, r);
U=U*sqrt(D);
V=V*sqrt(D);
P=[U;V];

W=eye(r);
[~,j] = find(A);
[rows,cols]=ind2sub([m,n], j);
tic;
for k=1:maxIter
    Pold=P;
    
    funU=@(x) objectiveU(x,V,W,A,b,beta);
    stepFunU=@(x,d) tu(x,V,W,A,b,beta,d);
    
    
    U=steepestDescent(U, funU, stepFunU, SDopts);
    
    funV=@(x) objectiveV(U,x,W,A,b,beta);
    stepFunV=@(x,d) tv(U,x,W,A,b,beta,d);
    
    V=steepestDescent(V, funV, stepFunV, SDopts);
    
    P=[U;V];
    
    [Vw,D]=eig(U'*U+V'*V);
    W=Vw*diag(f(diag(D),gamma))*Vw';
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
    if print_bool && k>1
        fprintf("it:  %i  rank:  %i  obj:  %f  rel: %f  gamma: %f \n", k, sum(diag(D)>0.0001), obj(k), abs(obj(k-1)-obj(k)), gamma)
    end

end

    function [f,g]=objectiveU(x,V,W,A,b,beta)
        f=trace(W*(x'*x))+(beta/2)*norm(A*vec(sparse_multiply(x,V, rows, cols,m,n))-b,'fro')^2;
        g=2*x*W+beta*reshape(A'*(A*vec(sparse_multiply(x,V, rows, cols,m,n))-b),[m,n])*V;
    end

    function [f,g]=objectiveV(U,x,W,A,b,beta)
        f=trace(W*(x'*x))+(beta/2)*norm(A*vec(sparse_multiply(U,x, rows, cols,m,n))-b,'fro')^2;
        g=2*x*W+ beta*reshape(A'*(A*vec(sparse_multiply(U,x, rows, cols,m,n))-b),[m,n])'*U;
    end


    function [t]=tu(U,V,W,A,b,beta,g)
        UV=sparse_multiply(U,V, rows, cols,m,n);
        gV=sparse_multiply(g,V, rows, cols,m,n);
        t=(2*trace(g'*U*W)+(beta)*(A*vec(gV))'*(A*(vec(UV))-b))/(beta*norm(A*vec(gV), 'fro')^2+2*trace(g'*g*W));
    end

    function t=tv(U,V,W,A,b,beta,g)
        
        
        UV=sparse_multiply(U,V, rows, cols,m,n);
        Ug=sparse_multiply(U,g, rows, cols,m,n);
        t=(2*trace(g'*V*W)+(beta)*(A*vec(Ug))'*(A*(vec(UV))-b))/(beta*norm(A*vec(Ug), 'fro')^2+2*trace(g'*g*W));
        
    end
end


