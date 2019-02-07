function [f,g]= testFunc(x,h)

f=trace(h*x*x');
g=2*h*x;

end