function X=calculateSpaseX(P,rows,cols)
[n,r]=size(P);
Pl=P(rows,:);
Pr=P(cols,:);
S=sum((Pl.*Pr), 2);
X=sparse(rows, cols, S, n,n);
end