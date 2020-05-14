load_path
clear
m=1000;
n=500;
r=15;
d=0.05;
p=0.2;


[D, A, b,M]=generateMatrixCompletionProblem(m,n,r,p,d);
Omega=(M~=0);
[~,j] = find(A);
[rows,cols]=ind2sub([m,n], j);
opts.r=2*r;
opts.maxIter=2000;
opts.fTol=1e-5;
opts.beta=0.0001;
gamma=norm(M,'fro')/(sqrt(sum(sum(Omega))/(m*n))*r);

opts.obj=@(U,V) norm(D-U*V','fro')/norm(D,'fro');
%opts.obj=@(U,V) sum(1- gamma./(gamma+eig([U;V]'*[U;V])))...
%        +(0.001/2)*norm(A*vec(sparse_multiply(U,V, rows, cols,m,n))-b)^2;

opts.ex=0;
opts.gamma=1;
[ U_sor,V_sor , obj_sor,time_sor] = genAltMin_v2(M,Omega,opts);
norm(D-U_sor*V_sor','fro')/norm(D,'fro')

Omega=(M~=0);
opts.beta=0.01;

[ U_ASD,V_ASD , obj_ASD,time_ASD] = GenASD(M,opts);
norm(D-U_ASD*V_ASD','fro')/norm(D,'fro')

opts.gamma=10;
opts.beta=0.01;
opts.maxIter=50;
[ U,V , obj,time] = GenAltMin(M,A,b,opts );
norm(D-U*V','fro')/norm(D,'fro')

figure()
loglog(time_ASD, obj_ASD,'-*', 'MarkerSize', 5,'LineWidth',1);
hold on
loglog(time_sor, obj_sor,'-o', 'MarkerSize', 5,'LineWidth',1);

loglog(time, obj,'-s','LineWidth',1);
xlabel("Time (s)")
ylabel("RFNE")
legend(["ASD","SOR", "AM"])
%legend(["ASD","AM"])
set(gca,'fontsize', 14);
