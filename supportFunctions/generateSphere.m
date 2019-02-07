function[x,y,z]=generateSphere(n)
%generate sphere in sherical coordinates
theta=pi*rand(n,1);
phi=2*pi*rand(n,1);
r=1;

%convert to cartisian

x=r*sin(theta).*cos(phi);
y=r*sin(theta).*sin(phi);
z=r*cos(theta);

%scatter3(x,y,z);
end