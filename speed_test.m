m=300;
n=200;

[D, A, b,M]=generateMatrixCompletionProblem(m,n,5,0.4, 0.1);
[~,j] = find(A);
[rows,cols]=ind2sub([m,n], j);



opts.r=10;
opts.maxIter=200;
opts.mu=0.01;
opts.xTol=1e-4;

gamma=norm(b)/sqrt(2*0.4*10);

opts.f=@(x) gamma./(gamma+x).^2;

opts.obj=@(U,V) sum(1- gamma./(gamma+eig([U;V]'*[U;V])))...
    +(opts.mu/2)*norm(A*vec(sparse_multiply(U,V, rows, cols,m,n))-b)^2;

opts.ex=0;

[ U_ASD,V_ASD , obj_ASD,time_ASD] = GenASD(M,A,b,opts );
norm(D-U_ASD*V_ASD','fro')/norm(D,'fro')
% opts.ex=0;
% [ U_ASD0,V_ASD0 , obj_ASD0,time_ASD0] = GenASD(M,A,b,opts );
%norm(D-U_ASD0*V_ASD0','fro')/norm(D,'fro')

opts.maxIter=200;
[ U,V , obj,time] = GenAltMin(M,A,b,opts );
norm(D-U*V','fro')/norm(D,'fro')

figure()
semilogy(time_ASD, obj_ASD,'-o', 'MarkerSize', 2,'LineWidth',1);
hold on
semilogy(time, obj,'-s','LineWidth',1);
xlabel("Time (s)")
ylabel("Objective")
legend(["ASD","AltMin"])
set(gca,'fontsize', 14);


% figure()
% loglog( obj_ASD);
% hold on
% loglog( obj_ASD0);
% xlabel("iteration")
% ylabel("Objective")
% legend(["ASD with Extrapolation Term","ASD without Extrapolation Term"])
