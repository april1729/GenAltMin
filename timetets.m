clear;
m=5000;n=1000;r=13;d=0.1;p=0.3;
%m=10;n=10;r=2;d=0;p=0.3;
[D, A, b,M]=generateMatrixCompletionProblem(m,n,r,p, d);
opts.r=2*r;
opts.mu=0.1;
opts.gamma=1e5;
opts.maxIter=500;
opts.f=@(x,gamma) 1000* gamma./(gamma+x).^2;
[rows,cols]=ind2sub([m,n], find(M~=0));
gamma=opts.gamma
opts.obj=@(U,V) sum(1- gamma./(gamma+eig([U;V]'*[U;V])))...
    +(opts.mu/2)*norm(A*vec(sparse_multiply(U,V, rows, cols,m,n))-b)^2;
Omega=(M~=0);
disp long
tic;[X1,Y1,Out] = lmafit_mc_adp(m,n,3*r,find(Omega),M(find(Omega)),[]);toc
norm(D-X1*Y1,'fro')/norm(D,'fro')


opts.mu=0.1;
% tic;[ U,V , obj,time] = GenASD(M,opts );toc
% norm(D-U*V','fro')/norm(D,'fro')

opts.maxIter=5000;

tic;[X,Y,obj,time] =  genAltMin_v2(M,Omega,opts);toc;
norm(D-X*Y','fro')/norm(D,'fro')

tic; [ U,V , obj,time] = GenASD(M,opts ); toc;
norm(D-U*V','fro')/norm(D,'fro')
