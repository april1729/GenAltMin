clear;
m=5000;n=1000;r=10;d=0.1;p=0.2;
[D, A, b,M]=generateMatrixCompletionProblem(m,n,r,p, d);
opts.r=2*r;
opts.mu=0.01;
opts.gamma=1000;
opts.maxIter=500;
opts.f=@(x,gamma)  gamma./(gamma+x).^2;
[rows,cols]=ind2sub([m,n], find(M~=0));

opts.obj=@(U,V) sum(1- opts.gamma./(opts.gamma+eig([U;V]'*[U;V])))...
    +(opts.mu/2)*norm(A*vec(sparse_multiply(U,V, rows, cols,m,n))-b)^2;

Omega=(M~=0);

% tic;[ U,V , obj,time] = GenASD(M,opts );toc
% norm(D-U*V','fro')/norm(D,'fro')

opts.maxIter=5000;
opts.scaled=1;
tic;[X,Y,obj,time] = GenASD(M,opts);toc
norm(D-X*Y','fro')/norm(D,'fro')

opts.scaled=0;
tic;[X0,Y0,obj0,time0] = GenASD(M,opts);toc
norm(D-X0*Y0','fro')/norm(D,'fro')


loglog(obj0);
hold on
loglog(obj);
legend(["ASD", "Scaled ASD"]);