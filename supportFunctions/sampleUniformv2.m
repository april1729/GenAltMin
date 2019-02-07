function [A, b, M] = sampleUniformv2(D,pavg)

[m,n]=size(D);
randMat=rand(m,n);
sampleIndex=find(randMat<=pavg);
[I,J] = ind2sub([m,n],sampleIndex);

A=sparse(1:length(I), sub2ind([m,n], I,J), ones(length(I), 1), length(I), n*m);
b=D(sub2ind([m,n], I,J));

M=sparse(I,J, b, m,n);
end