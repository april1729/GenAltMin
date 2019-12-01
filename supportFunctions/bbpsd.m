function [ P , obj] = bbpsd(P0,W,A,b,params )
%BBALPSD Summary of this function goes here
%   Detailed explanation goes here
%   This solves the problem:
%
%   min trace(W*P'*P)+mu/2 ||A*vec(P*P')-b||^2

[n,r]=size(P0);
mu=params{3};
maxIter=params{4};
eps=params{5};
[~,j] = find(A);
[rows,cols]=ind2sub([n,n], j);
P=P0;
X=calculateSpaseX(P,rows,cols);
s=2*P*W+2*mu*reshape(A'*(A*vec(X)-b),[n,n])*P;
Xv=vec(X);
AtA=A'*A;
Atb=A'*b;

f=@(x) trace(W*x'*x)+mu/2*(norm(A*vec(x*x')-b))^2;

for k=1:maxIter
    sold=s;
    s=2*P*W+2*mu*reshape((AtA*Xv-Atb),[n,n])*P;
    if k==1
        alpha=fminsearch(@(a) f(P-a*s), 1e-3);
    else
        alpha=((vec(P-Pold)'*vec(s-sold))/(vec(s-sold)'*vec(s-sold)));
        %alpha=10e-5;
        if alpha<0
            alpha=abs(alpha)/10;
        end
%        alpha=fminsearch(@(a) f(P-a*s), 1e-3);
        
    end
    Pold=P;
    
    P=P-alpha*s;
    X=calculateSpaseX(P,rows,cols);
    Xv=vec(X);
    
    
    obj(k)=trace(W*P'*P)+(mu/2)*norm(A*Xv-b,'fro')^2;
    
    if k>2
        if abs(log(obj(k-1))-log(obj(k-2)))<eps%norm((P-Pold), 'fro')/norm(Pold, 'fro')<eps
            break
        end
    end
    
end
X=P*P';

end
