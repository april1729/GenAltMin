function [ X, S ] = sparsityRegulatedRankMinimizationConvex( X, A,b, delta, E,g, params)
%SPARSITYREGULATEDRANKMINIMIZATIONCONVEX Summary of this function goes here
%   Detailed explanation goes here
[m,n]=size(X);
gamma=params{1};
cvx_begin quiet
variable X(m,n)
variable S(m,n)
minimize norm_nuc(X)+gamma*sum(sum(abs(S)))
subject to
if ~isempty(A)
    square_pos(norm(A*vec(X+S)-b)) <= delta
end
if ~isempty(E)
    E*vec(X+S)==g;
end
cvx_end

end

