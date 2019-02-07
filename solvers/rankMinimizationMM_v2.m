function [X,obj]=rankMinimizationMM_v2(x0, A, b, delta, E, g, params)
% Params start with a string, either MR, SCAD, or LogDet
[m,n]=size(x0);
eps=0.001;
maxIter=1000;
if isempty(params)
    f=@(x) 1-(0.1./(0.1+x));
    df=@(x) 0.1./(0.1+x).^2;
else
    f=params{1};
    df=params{2};
end

X=x0;

for k=1:maxIter
    Xold=X;
    [U,S,V]=svd(X,'econ');
    W=U*diag(df(diag(S)))*V';
try    
    cvx_begin quiet
    cvx_precision low
    variable X(m,n)
    minimize trace(transpose(X)*W))
    subject to
    if ~isempty(A)
        square_pos(norm(A*vec(X)-b)) <= delta
    end
    if ~isempty(E)
        E*vec(X)==g;
    end
    cvx_end
catch
    disp("oh no")
end
    obj(k)=sum(f(svd(X)));    
    if norm(X-Xold,'fro')/norm(X,'fro')<eps
        break;
    end
end