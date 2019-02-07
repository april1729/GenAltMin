percents=[0.125,0.25,0.5,1,2]./100;
methods={@(x) 1,@(x) scad(x, 50,2)/scad(0,50,2),@(x) 50^2./(50+x).^2, @(x) 50./(50+x),  @(x) (10+x).^(-0.5)/(10^-0.5)};
[pt, F]=ReadOFF("sphere1k.off");
%load("cow.mat")
d=3;
r=10;
Dist=calcDistance(pt*pt');
n=length(Dist);
P0=rand(n,r);
eigsCellArray={};
opts.exact=false;
noise=randn(n,n);
noise=0.5*(noise+noise');
opts.mu=1;
for p=1:length(percents)
    [A, b] = sampleUniformGram(Dist+0.05*noise,percents(p));
    for m=1:length(methods)
        opts.f=methods{m};
        fprintf("Percent: %i method: %i \n", p, m)
        tic;[X,obj]=gnrtn(P0,A,b,opts);toc;
        time=toc;
        ptr=reconstructPoints(X, eye(3));
        figure()
        ViewMesh(ptr,F);
        filenameFig=char("spheresNoise/m"+m+"p"+p+".fig");
        filenamePNG=char("spheresNoise/m"+m+"p"+p+".png");
        saveas(gcf,filenameFig)
        saveas(gcf,filenamePNG)
        title("method "+ m+" with "+100*percents(p)+ " percent of data")
        errorsArray(m,p)=calcErrorPoints(X, pt);
        timeArray(m,p)=time;
        nextEigArray(m,p)=min(eigs(X,4));
        rankArray(m,p)=sum(eigs(X,10)>0.0001);
        eigsCellArray{m,p}=eigs(X,10);
    end
end
save("spheresNoise/spheresNoise.mat")


