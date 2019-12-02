function [X,Y,obj,time] = genAltMin_v2(M,Omega,opts)
if isfield(opts, 'mu')
    mu=opts.mu;
else
    mu=1;
end
if isfield(opts, 'f')
    f=opts.f;
else
    f=@(x) 50^2./(50+x).^2;
end
if isfield(opts, 'maxIter')
    maxIter=opts.maxIter;
else
    maxIter=100;
end
if isfield(opts, 'xTol')
    xTol=opts.xTol;
else
    xTol=1e-4;
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
[m,n]=size(M);
[rows,cols]=ind2sub([m,n], find(Omega==1));
Known=find(Omega==1);
data=M(Known);
% [U, D, V]=svds(M, r);
% X=U*sqrt(D);
% Y=V*sqrt(D);
X=randn(m,r);
Y=randn(n,r);
P=[X;Y];
theta=1;
W=eye(r);
Z=M;
tic;
time=[];obj=[];
for k=1:maxIter
    
    Xp=mu*(Z*Y)*pinv((W+mu*Y'*Y));
    X=theta*Xp+(1-theta)*X;
    
    Yp=mu*(Z'*X)*(pinv(W+mu*X'*X));
    Y=theta*Yp+(1-theta)*Y;
    Z = X*Y';  Res =  data - Z(Known);
    Z(Known) = data  + theta * Res;
    [Vw,D]=eig(X'*X+Y'*Y);
    
    
    W=Vw*diag(f(diag(D)))*Vw';
    
    Pold=P;
    P=[X;Y];
    
    
    
    time(k)=toc;
        if isfield(opts, "obj")
            obj(k)=sum(diag(D)>0.1)+mu*norm(Res)^2;
            if k>2
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
