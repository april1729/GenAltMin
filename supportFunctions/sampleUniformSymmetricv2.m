function [A, b, M] = sampleUniformSymmetricv2(D,pavg)
[n,~]=size(D);
randMat=rand(n,n);
sampleIndex=find(randMat<=pavg);
[I,J] = ind2sub([n,n],sampleIndex);
I=I';
J=J';
isUpper=I<J;
isDiag=(I==J);
Iold=I;
I=[I(isUpper),J(isUpper), I(isDiag)];
J=[J(isUpper),Iold(isUpper), J(isDiag)];

A=sparse(1:length(I), sub2ind([n,n], I,J), ones(length(I), 1), length(I), n^2);
b=D(sub2ind([n,n], I,J))';

M=sparse(I,J, b, n,n);
end