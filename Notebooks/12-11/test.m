
X=rand(10,5);

f=@(x) 1./(1+x);
g=@(x) -1./(1+x).^2;
F=@(X) sum(1./(1+eig(X*X')));


% Numerical gradient
[V,D]=eig(X*X');
df=V*diag(g(eig(X*X')))*V';

%analytic gradient
dx=0.001;
for i =1:10
    for j=1:5
        dX=zeros(10,10);
        dX(i,j)=dx;
        grad(i,j)=(F(X+dX)-F(X))/dx;
    end
end

norm(df-grad)