function [X,obj]=rankMinimizationPSDMM(x0, A, b, delta, E, g, params)
% Params start with a string, either MR, SCAD, or LogDet
n=length(x0);

eps=0.001;
maxIter=1000;

if isempty(params)
    method="MR";
    lambda=0.1;
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
        W=inv((1/lambda)*X+eye(n))^2;
    elseif strcmp(method, "LogDet")
        W=inv((1/lambda)*X+eye(n))^1;
    elseif strcmp(method, "SCAD")
        [V,D]=eig(X);
        [~,w]=scad(diag(D),beta,theta);
        W=real(V*diag(w)*V');
    end
        
    
    Xold=X;
try    
    cvx_begin quiet
    cvx_precision low
    variable X(n,n)
    minimize abs(sum(sum(W.*X)))
    subject to
    if ~isempty(A)
        square_pos(norm(A*vec(X)-b)) <= delta
    end
    if ~isempty(E)
        E*vec(X)==g;
    end
    X == semidefinite(n)
    cvx_end
catch
    disp("oh no")
end
    X=full(X);
    x=eig(X);
    if strcmp(method, "MR")
        obj(k)=n-sum(lambda./(x+lambda));
    elseif strcmp(method, "LogDet")
        obj(k)=sum(log(x/lambda+1));
    elseif strcmp(method, "SCAD")
        obj(k)=sum(scad(x, beta, theta));
    end
    

    
    if norm(X-Xold,'fro')/norm(X,'fro')<eps
        break;
    end
    
end