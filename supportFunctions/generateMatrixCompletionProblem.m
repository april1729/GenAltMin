function [D, A, b,M]=generateMatrixCompletionProblem(m,n,r,p, noiseLevel)
D=randn(m,r)*randn(r,n);
noise=randn(m,n);
D_noise=D+noiseLevel*noise;
[A,b,M] = sampleUniform(D_noise,p);
end