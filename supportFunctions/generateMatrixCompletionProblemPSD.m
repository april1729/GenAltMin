function [D, A, b,M]=generateMatrixCompletionProblemPSD(n,r,p, noiseLevel)
x=randn(r,n);
noise=randn(n);
noise=noiseLevel*(noise'+noise);

D=x'*x;
[A,b,M] = sampleUniformSymmetricv2(x'*x+noise,p);
end