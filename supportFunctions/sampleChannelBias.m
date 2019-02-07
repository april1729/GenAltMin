function [M, b ,row, col ] = sampleChannelBias( A, channels, pavg, q)
%CORRELATEDSAMPLING Summary of this function goes here
%   Detailed explanation goes here

[m,n]=size(A);

for i=1:m 
   for j=1:length(channels)
       if rand()>q
           if j==1
               A(i, 1:channels(j))=0;
           else
               A(i, (channels(j-1)+1):channels(j))=0;
           end
       end
   end
end

M=A.*(rand(m,n)>pavg);
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