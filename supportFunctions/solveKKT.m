function X=solveKKT(W, M)

[m,n]=size(M);
X=M;
P=(M~=0);

for i =1:m
    
    S= find(P(i,:)==0);
    Sc= find((1-P(i,:))==0);

    b=-W(:, Sc)*M(i,Sc)';
    A=W(:, S);
    
    %x=pinv(A)*b;
    x=A\b;
    X(i, S)=x;
    
end