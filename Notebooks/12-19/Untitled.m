%% oh no this isnt actually concave but its fine i can fix it

[D, A, b]=generateMatrixCompletionProblem(100,50,3,0.5, 0);

lambda=1;
f=@(x) 1-(lambda./(sqrt(x+lambda^2)));
g=@(x) (1/2)*lambda./(sqrt(x+lambda^2)).^(3/2);

[X,obj]=rankMinimizationMM(zeros(100,50), [],[],0, A,b, {"MR", lambda});
[Xsqrt,objsqrt]=rankMinimizationMM(zeros(100,50), [],[],0, A,b, {f,g});

plot(svd(X))
hold on
plot(svd(Xsqrt))

norm(X-D)/norm(D)
norm(Xsqrt-D)/norm(D)

%% so that works and all, but it still isnt concave!!! 
% What if we just do this:  let 
%
% $$f(x)=1-\frac{\lambda}{\sqrt{x}+\lambda} $$ 
%
% and then we have 
%
% $$f'(x)=\frac{1}{2} x^\frac{-1/2} \frac{\lambda}{(\sqrt{x}+\lambda)^2} $$
%
% But, this is undeefined at x=0, so we just use min(f(x), 10000) instead.

lambda=.5;  

f=@(x) 1-lambda./(sqrt(x+0.01)+lambda);
g=@(x) 0.5*(x+0.01).^(-0.5).*(lambda./(sqrt(x+0.01)+lambda));

[D, A, b]=generateMatrixCompletionProblem(100,50,3,0.5, 0);


[X,obj]=rankMinimizationMM(zeros(100,50), [],[],0, A,b, {"MR", lambda});
[Xsqrt,objsqrt]=rankMinimizationMM(zeros(100,50), [],[],0, A,b, {f,g});
figure()
plot(svd(X))
hold on
plot(svd(Xsqrt))

norm(X-D)/norm(D)
norm(Xsqrt-D)/norm(D)

%% we can't just shift it over

% What if we just do this:  let 
%
% $$f(x)=1-\frac{\lambda}{\sqrt{x+\epsilon}+\lambda} $$ 
%
% but its always going to be a function of $x^2$, so in the end it looks
% like 
% $$f(x^2)=1-\frac{\lambda}{\sqrt{x^2+\epsilon}+\lambda} $$ 
%
% As opposed to 
% $$f(x^2)=1-\frac{\lambda}{\sqrt{x+\epsilon}^2+\lambda} $$ 
% but lets try it non the less.  

f=@(x) 1-lambda./(sqrt(x+0.01)+lambda);
g=@(x) 0.5*(x+0.01).^(-0.5).*(lambda./(sqrt(x+0.01)+lambda));

[D, A, b]=generateMatrixCompletionProblem(100,50,3,0.5, 0);


[X,obj]=rankMinimizationMM(zeros(100,50), [],[],0, A,b, {"MR", lambda});
[Xsqrt,objsqrt]=rankMinimizationMM(zeros(100,50), [],[],0, A,b, {f,g});
figure()
plot(svd(X))
hold on
plot(svd(Xsqrt))

norm(X-D)/norm(D)
norm(Xsqrt-D)/norm(D)
