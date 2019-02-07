function mu=coherence(M)

[m,n]=size(M);
if m<n
    M=M';
    [m,~]=size(M);
end


[u,~,~]=svds(M,rank(M));
vecnorm=zeros(m,1);
for i=1:m
vecnorm(i)=norm(u(i,:))^2;
end
mu=max(vecnorm);