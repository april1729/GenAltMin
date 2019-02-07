%% Matrix Completion for Multi-Label Image Classification
% DESCRIBE IDEA HERE


%% Load the data
% In a seperate script i tool the MSRC Multi label image dataset and gor it
% into the form described in the paper.  

load("MSRC_Histograms_Matrix")
[m,n]=size(Z);
%% Modeling Strategy 1: Just hold it all equal
testImages=datasample(1:n, 80);
numSampled=0;
numNotSampled=0;
Aeq=sparse(m*n,m*n);
A=sparse(80*14,m*n);

for i=1:m
    for j=1:n
        if (~ismember(j, testImages) || ~ismember(i, 1:14))
            numSampled=numSampled+1;
            Aeq(numSampled, (j-1)*m+i)=1;
            beq(numSampled,1)=Z(i,j);
        else
            numNotSampled=numNotSampled+1;
            A(numNotSampled, (j-1)*m+i)=1;
        end
    end
end
Aeq=Aeq(1:numSampled,:);
norm(Aeq*vec(Z)-beq)


%% Does it actually work?
% Lets find out.

eps=0.01;
maxIter=1000;
lambda=100;
X=zeros(m,n);

for k=1:maxIter
    W=inv((1/lambda)*(X'*X)+eye(n))^2;
    Xold=X;

    cvx_begin quiet
    cvx_precision low
    variable X(m,n)
    minimize norm(X*W^0.5, 'fro')
    subject to
    Aeq*vec(X)==beq;    
    A*vec(X)>=0;
    A*vec(X)<=1;
    cvx_end

    X=full(X);
    obj(k)=n-trace(inv(X'*X/lambda+eye(n)));
    if norm(X-Xold,'fro')/norm(X,'fro')<eps
        break;
    end
    
end


ytest=Z(1:14,testImages);
ypredict=(X(1:14,testImages)>0.9);

sum(sum(ypredict==ytest))/(14*80)