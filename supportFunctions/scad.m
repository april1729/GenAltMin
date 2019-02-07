function [g,f]=scad(x, beta, theta)
f=zeros(size(x));
g=zeros(size(x));

for i =1:length(x)
    if x(i)<= beta
        f(i)=beta*x(i);
        g(i)=beta;
    elseif x(i)<= beta*theta
        f(i)=(-x(i)^2+2*beta*theta*x(i)-beta^2)/(2*(theta-1));
        g(i)=(-2*x(i)+2*beta*theta)/(2*(theta-1));
    else
        f(i)=beta^2*(1+theta)/2;
        g(i)=0;
    end
end

end
