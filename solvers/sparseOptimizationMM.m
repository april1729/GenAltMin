function [x,obj]=sparseOptimizationMM(x, A, b, params)
% Params start with a string, either MR, SCAD, or LogDet
n=length(x);
eps=0.00001;
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

for k=1:maxIter
    if strcmp(method, "MR")
        w=lambda./(x+lambda).^2;
    elseif strcmp(method, "LogDet")
        w=((1/lambda)*x+1).^-1;
    elseif strcmp(method, "SCAD")
        [~,w]=scad(x, beta, theta);
    end
    
    xold=x;
    x = linprog(w,[],[],A,b, zeros(n,1), []);
    
    if strcmp(method, "MR")
        obj(k)=n-sum(lambda./(x+lambda));
    elseif strcmp(method, "LogDet")
        obj(k)=sum(log(x/lambda+1));
    elseif strcmp(method, "SCAD")
        obj(k)=sum(scad(x, beta, theta));
    end
    
    if norm(x-xold,'fro')/norm(x,'fro')<eps
        break;
    end
end

end