function [D, A, b]=generateMatrixCompletionProblem(m,n,r,p, noiseLevel)
D=randn(m,r)*randn(r,n);
noise=randn(m,n);
D=D+noiseLevel*noise;
[A,b] = sampleUniform(D,p);
end