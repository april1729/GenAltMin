function X=steepestDescent(X0, objFun, stepFun, opts)
% Minimizes an unconstrained optimization problem by steepest descent.
%
% X0 - starting pint
% objFun - Function that takes one arguement, the current iterate and returns the objective and gradient
%			E.X. [f,g]=objFun(X);
% setFun - Function that takes the current iterate and gradient at that iterate and returns the optimal steep size
%			e.x. [t]=stepFun(X,g)

if isfield(opts, 'maxIter')
    maxIter=opts.maxIter;
else
    maxIter=100;
end
if isfield(opts, 'fTol')
    fTol=opts.fTol;
else
    fTol=1e-4;
end
if isfield(opts, 'xTol')
    xTol=opts.xTol;
else
    xTol=1e-4;
end
X=X0;
[f,~]=objFun(X);


for i=1:maxIter
    Xold=X;
    fOld=f;
    [f,g]=objFun(X);
    X=X-stepFun(X,g)*g;
    obj(i)=f;
    if norm(X-Xold, 'fro')/norm(Xold, 'fro')<xTol && abs(f-fOld)/abs(fOld) <fTol
        break
    end
end

end
