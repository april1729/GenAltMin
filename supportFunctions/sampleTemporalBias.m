function [M, b ,row, col ] = sampleTemporalBias( A, pavg, tau)
%TEMPORALCORRELATION Summary of this function goes here
%   Detailed explanation goes here

[m,n]=size(A);
M=A;
for j=1:n
    for i=1:m
        if rand()<pavg
            M(i:(min(m, i+tau-1)), j)=0;
        end
    end
end



b=[];
row=[];
col=[];
for i =1:m
    for j=1:n
        if M(i,j)~=0
            b=[b, M(i,j)];
            row=[row, i];
            col=[col, j];
        end
    end
end

end

