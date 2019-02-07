percents=[0.5, 1, 2.5,5,10]./100;
methods={@(x) 1, @(x) 10./(10+x), @(x) 10./(10+x).^2, @(x) (10+x).^(-0.5),@(x) scad(x, 5, 2)};
m=1000;
n=500;
r=10;
eigsCellArray={};
for p=1:length(percents)
    [D, A, b]=generateMatrixCompletionProblem(m,n,r,percents(p), 0);
    
    for meth=1:length(methods)
        [X,obj]=rankMinimizationMM_v2(zeros(m,n), [],[],0, A,b, {@(x) 1,methods{meth}});
        time=toc;
        ptr=reconstructPoints(X, eye(3));
        figure()
        ViewMesh(ptr,F);
        title("method "+ meth+" with "+100*percents(p)+ " percent of data")
        errorsArray(meth,p)=calcErrorPoints(X, pt);
        timeArray(meth,p)=time;
        nextEigArray(meth,p)=min(eigs(X,4));
        rankArray(meth,p)=sum(eigs(X,10)>0.5);
        eigsCellArray{m,p}=eigs(X,10);
    end
end