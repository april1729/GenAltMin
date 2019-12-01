function [X]=sparse_multiply(U,V, rows, cols,m,n)

Us=U(rows,:);
Vs=V(cols,:);
S=sum(Us.*Vs, 2);
X=sparse(rows, cols, S, m,n);



