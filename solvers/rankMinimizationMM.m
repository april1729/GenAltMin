function [X,obj]=rankMinimizationMM(x0, A, b, delta, E, g, params)
% Params start with a string, either MR, SCAD, or LogDet
[m,n]=size(x0);
eps=0.01;
maxIter=1000;

if isempty(params)
    method="MR";
    lambda=0.1;
    
elseif ~(isstring(params{1}) || ischar(params{1}))
    f=params{1};
    df=params{2};
    method="other"
elseif strcmp("MR", params{1})
    lambda=params{2};
    method="MR";
elseif strcmp("LogDet", params{1})
    lambda=params{2};
    method="LogDet";
elseif strcmp("SCAD", params{1})
    method="SCAD";
    beta=params{2};
    theta=params{3};
end

X=x0;

for k=1:maxIter
    if strcmp(method, "MR")
        W=inv((1/lambda)*(X'*X)+eye(n))^2;
    elseif strcmp(method, "LogDet")
        W=inv((1/lambda)*(X'*X)+eye(n))^1;
    elseif strcmp(method, "SCAD")
        [P,D]=eig((X'*X));
        [~,w]=scad(diag(D),beta,theta);
        W=real(P*diag(w)*P');
    elseif strcmp(method, "other")
        [P,D]=eig((X'*X));
        w=df(diag(D));
        W=real(P*diag(w)*P');
    end
        
    
    Xold=X;
try    
    cvx_begin quiet
    cvx_precision low
    variable X(m,n)
    minimize norm(X*W^0.5, 'fro')
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
    X=full(X);
    x=eig((X'*X));
    if strcmp(method, "MR")
        obj(k)=n-sum(lambda./(x+lambda));
    elseif strcmp(method, "LogDet")
        obj(k)=sum(log(x/lambda+1));
    elseif strcmp(method, "SCAD")
        obj(k)=sum(scad(x, beta, theta));
    elseif strcmp(method, "other")
        obj(k)=sum(f(x));
    end
    

    
    if norm(X-Xold,'fro')/norm(X,'fro')<eps
        break;
    end
    
end