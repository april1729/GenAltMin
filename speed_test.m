clear
m=1000;
n=500;
r=10;
d=0.3;
p=0.1;


[D, A, b,M]=generateMatrixCompletionProblem(m,n,r,p,d);
Omega=(M~=0);
[~,j] = find(A);
[rows,cols]=ind2sub([m,n], j);
opts.r=2*r;
opts.maxIter=400;
gamma=norm(M,'fro')/(sqrt(sum(sum(Omega))/(m*n))*r);

%opts.obj=@(U,V) norm(D-U*V','fro')/norm(D,'fro');
opts.obj=@(U,V) sum(1- gamma./(gamma+eig([U;V]'*[U;V])))...
        +(0.001/2)*norm(A*vec(sparse_multiply(U,V, rows, cols,m,n))-b)^2;



[ U_sor,V_sor , obj_sor,time_sor] = genAltMin_v2(M,Omega,opts);
norm(D-U_sor*V_sor','fro')/norm(D,'fro')

Omega=(M~=0);
opts.beta=0.001;

[ U_ASD,V_ASD , obj_ASD,time_ASD] = GenASD(M,opts);
norm(D-U_ASD*V_ASD','fro')/norm(D,'fro')


opts.gamma=100;

opts.maxIter=500;
%[ U,V , obj,time] = GenAltMin(M,A,b,opts );
%norm(D-U*V','fro')/norm(D,'fro')

figure()
semilogy(time_ASD, obj_ASD,'-o', 'MarkerSize', 2,'LineWidth',1);
hold on
semilogy(time_sor, obj_sor,'-o', 'MarkerSize', 2,'LineWidth',1);

%semilogy(time, obj,'-s','LineWidth',1);
xlabel("Time (s)")
ylabel("RFNE")
legend(["Alternating Steepest Descent","SOR", "Alternating Minimization"])
set(gca,'fontsize', 14);
