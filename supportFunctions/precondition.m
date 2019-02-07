function [ M ] = preconditionMatrix( M )
%PRECONDITION_MATRIX Summary of this function goes here
%   Detailed explanation goes here
    M=M-mean((M));
    M=M./sqrt(var(M));
end

