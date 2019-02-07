function M=ICMC(row, col, b,m,n,r)

G=graph();
G=addnode(G, m+n);
for i=1:length(row)
    G=addedge(G, row(i), m+col(i));
end

X=zeros(m,r);
Y=zeros(n,r);



for i = 1:m
    if degree(G,i+m)>r
        infectedU=neighbors(G,i+m)
        infectedU=infectedU(1:r)'
        break
    end
end
    
infectedV=[];
X(infectedU,:)=eye(r,r);


oldNumInfected=0;
while (length(infectedU)+length(infectedV)<m+n) && (oldNumInfected<(length(infectedU)+length(infectedV)))
    oldNumInfected=length(infectedU)+length(infectedV)
    for vj = setdiff(1:n, infectedV)
        s=intersect(neighbors(G, vj+m), infectedU);
        if length(s)>=r
            Y(vj, :)=pinv(X(s,:))*b(s)';
            infectedV=[infectedV, vj];
            %break
        end
    end
    for uj = setdiff(1:m, infectedU)
        s=intersect(neighbors(G, uj)-m, infectedV);
        if length(s)>=r
            X(uj, :)=pinv(Y(s,:))*b(s)';
            infectedU=[infectedU, uj];
            %break
        end
    end
end

M=X*Y';
