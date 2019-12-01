function [D,S, M,A, b]=generateRobustMatrixCompletionProblem(m,n,r,s,p, noiseLevel)
D=randn(m,r)*randn(r,n);
S=zeros(m,n);
for i=datasample(1:n*m,s,'Replace',false)
    S(i)=randn();
    S(i)=S(i)+0.5*sign(S(i));
end
noise=randn(m,n);
M=D+S+noiseLevel*noise;
[A,b] = sampleUniform(M,p);
end