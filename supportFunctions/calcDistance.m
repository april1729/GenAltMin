function [ D ] = calcDistance( B )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
n=length(B);
        for i=1:n
            for j=1:n
                D(i,j)=B(i,i)+B(j,j)-2*B(i,j);
            end
        end



end

