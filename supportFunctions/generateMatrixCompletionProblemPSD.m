function [D, A, b]=generateMatrixCompletionProblemPSD(n,r,p, noiseLevel)
x=randn(r,n);
%noise=randn(n);
%noise=noiseLevel*(noise'+noise);

D=x'*x;
[A,b] = sampleUniformSymmetricv2(x'*x,p);
end