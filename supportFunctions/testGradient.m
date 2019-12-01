function testGradient(x,fun)
[n,r]=size(x);
grad=zeros(n,r);
dx=(sqrt(eps));
[f,g]=fun(x);
gtest=zeros(n,r);
for ind=datasample(1:n*r, round(log2(n*r)), 'Replace', false)
    [i,j]=ind2sub([n,r],ind);
    xt=x;
    xt(i,j)=xt(i,j)+dx;
    [ftemp,~]=fun(xt);
    grad(i,j)=(ftemp-f)/(dx);
    gtest(i,j)=g(i,j);
end

fprintf("Average error in gradient: %f \n",norm(gtest-grad,'fro')/round(log2(n*r)))

