function [points, distances, A, b]=generateSensorProblem(n,d,p, noise)
points=rand(n,d);
B=points*points';
distances=calcDistance(B);

b=[];
l=round(p*n*(n+1)/2);
isSampled=zeros(size(distances));
M=zeros(size(distances));
numSampled=0;
while numSampled<l
    i=randi([1,n]);
    j=randi([1,n]);
    while (isSampled(i,j)||i==j)
        i=randi([1,n]);
        j=randi([1,n]);
    end
    
    numSampled=numSampled+1;
    
    isSampled(i,j)=1;
    isSampled(j,i)=1;
    
    temp=zeros(n, n);
    temp(i, j)=-1;
    temp(j,i)=-1;
    temp(i, i)=1;
    temp(j,j)=1;
    A(numSampled, :)=vec(temp)';
    b(numSampled)=distances(i,j);
end
b=b';
end
