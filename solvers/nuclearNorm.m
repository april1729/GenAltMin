function [ X ] = nuclearNorm(x0,A,b,delta,E,g, params)
%NUCLEARNORMPSD Summary of this function goes here
%   Detailed explanation goes here

[m,n]=size(x0);



cvx_begin quiet
variable X(m,n)
minimize norm_nuc(X)
subject to  
if ~isempty(A)
square_pos(norm(A*vec(X)-b)) <= delta
end
if ~isempty(E)
E*vec(X)==g;
end
cvx_end


end

