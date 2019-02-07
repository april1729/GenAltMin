function [A, b, M] = sampleUniformSymmetric(D,pavg)
    [m,n]=size(D);
    l=round(pavg*n*(n+1)/2);
    datasample(1:10,5, 'Replace', false)
    
    isSampled=zeros(size(D));
    M=zeros(size(D));
    numSampled=0;
    while (sum(sum(M~=0)))<l
        i=randi([1,m]);
        j=randi([1,n]);
        while (isSampled(i,j)||i==j)
            i=randi([1,m]);
            j=randi([1,n]);
        end
        
        numSampled=numSampled+1;
        
        isSampled(i,j)=1;
        isSampled(j,i)=1;

        M(i,j)=D(i,j);
        temp=zeros(m, n);
        temp(i,j)=1;
        A(numSampled, :)=vec(temp)';
        b(numSampled)=D(i,j);
    end
    b=b';
end