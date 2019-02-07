function [ X, S, r_obj, s_obj ] = sparsityRegulatedRankMinimizationMM( X, A,b, delta, E,g, params)
%SPARSITYREGULATEDRANKMINIMIZATIONMM Summary of this function goes here
%   Detailed explanation goes here
[m,n]=size(X);
S=zeros(m,n);
gamma=params{1};
s_params=params{2};
r_params=params{3};
if isempty(s_params)
    s_method="lpcc";
    s_lambda=0.1;
elseif strcmp("lpcc", s_params{1})
    s_lambda=s_params{2};
    s_method="lpcc";
elseif strcmp("LogDet", s_params{1})
    s_lambda=s_params{2};
    s_method="LogDet";
elseif strcmp("SCAD", s_params{1})
    s_method="SCAD";
    s_beta=s_params{2};
    s_theta=s_params{3};
end


if isempty(r_params)
    r_method="MR";
    r_lambda=0.1;
elseif strcmp("MR", r_params{1})
    r_lambda=r_params{2};
    r_method="MR";
elseif strcmp("LogDet", r_params{1})
    r_lambda=r_params{2};
    r_method="LogDet";
elseif strcmp("SCAD", r_params{1})
    r_method="SCAD";
    r_beta=r_params{2};
    r_theta=r_params{3};
end


eps=0.01;
maxIter=50;

for k=1:maxIter
    if strcmp(s_method, "lpcc")
        s_W=s_lambda./(vec(abs(S))+s_lambda).^2;
    elseif strcmp(s_method, "LogDet")
        s_W=((1/s_lambda)*vec(abs(S))+1).^-1;
    elseif strcmp(s_method, "SCAD")
        [~,s_W]=scad(vec(abs(S)), s_beta, s_theta);
    end
    
    if strcmp(r_method, "MR")
        r_W=inv((1/r_lambda)*(X'*X)+eye(n))^2;
    elseif strcmp(r_method, "LogDet")
        r_W=inv((1/r_lambda)*(X'*X)+eye(n))^1;
    elseif strcmp(r_method, "SCAD")
        [P,D]=eig((X'*X));
        [~,d]=scad(diag(D),r_beta,r_theta);
        r_W=real(P*diag(d)*P' );
    end
    
    Xold=X;
    Sold=S;
try    
    cvx_begin quiet
    cvx_precision low
    variable X(m,n)
    variable S(m,n)
    minimize square_pos(norm(X*r_W^0.5, 'fro'))+gamma*transpose(s_W)*vec(abs(S))
    subject to
    if ~isempty(A)
        square_pos(norm(A*vec(X+S)-b)) <= delta
    end
    if ~isempty(E)
        E*vec(X+S)==g;
    end
    cvx_end
catch
    disp("oh no")
end
    X=full(X);
    S=full(S);
    
    
    x=eig((X'*X));
    if strcmp(r_method, "MR")
        r_obj(k)=n-sum(r_lambda./(x+r_lambda));
    elseif strcmp(r_method, "LogDet")
        r_obj(k)=sum(log(x/r_lambda+1));
    elseif strcmp(r_method, "SCAD")
        r_obj(k)=sum(scad(x, r_beta, r_theta));
    end
    
    if strcmp(s_method, "lpcc")
        s_obj(k)=sum(1-s_lambda./(vec(abs(S))+s_lambda));
    elseif strcmp(s_method, "LogDet")
        s_obj(k)=sum(log(vec(abs(S))/s_lambda+1));
    elseif strcmp(s_method, "SCAD")
        s_obj(k)=sum(scad(vec(abs(S)), s_beta, s_theta));
    end
    

    
    

    
   if norm(X-Xold,'fro')/norm(X,'fro')<eps
        break;
    end
    
end