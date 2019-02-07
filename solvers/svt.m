function X=svt(row, col, b, p, m, tau, delta, eps, maxIter)
Y=zeros(p,m);
for l=1:maxIter
    [U,S,V]=svd(Y);
    S=max(0, S-tau*eye(size(S)));
    X=U*S*V';
    Y=Y+delta*proj(X, row, col, b);
    if norm(proj(X, row, col, b))/norm(b)<eps
        break
    end
    %fprintf("iteration %i \t error %f \n", l, norm(proj(X, row, col, b))/norm(b))
end
end
function M=proj(X, row, col, b)
M=zeros(size(X));
for i=1:length(row)
    M(row(i), col(i))=b(i)-X(row(i), col(i));
end
end