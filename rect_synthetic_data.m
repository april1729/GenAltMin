load_path
clear
n=1000;
m=500;
r=10;
noiseLevel=.1;
numRuns=1;

p_range=0.2:0.1:1;

opts.r=2*r;
opts.maxIter=2000;
opts.xTol=1e-7;
opts.beta=1e-3;
gamma=0.1;
opts.gamma=gamma;

for run=1:numRuns
    for i=1:length(p_range)
        [D, A, b,M]=generateMatrixCompletionProblem(m,n,r,p_range(i), noiseLevel);
        Omega=(M~=0);
        [~,j] = find(A);
        [rows,cols]=ind2sub([m,n], j);
        opts.obj=@(U,V) sum(1- gamma./(gamma+eig([U;V]'*[U;V])))...
            +(0.001/2)*norm(A*vec(sparse_multiply(U,V, rows, cols,m,n))-b)^2;
        
        [ U_ti,V_ti , obj_nuc] = GenASD(M,opts );
        error(i, 1,run)=norm(U_ti*V_ti'-D,'fro')/norm(D, 'fro');
        rank_list(i,1,run)=sum(svd(U_ti*V_ti')>0.0001);
        
        [X]=svt(M, Omega,1, 2*r);
        
        error(i, 2, run)=norm(X-D,'fro')/norm(D, 'fro');
        rank_list(i,2,run)=sum(svd(X)>0.0001);
        
    end
end

figure()
plot(100*p_range, mean(error,3),'linewidth', 4)
hold on
plot(100*p_range, ones(size(p_range))*norm(A*vec(D)-b)/norm(b),'--','linewidth', 4)
xlabel("Percentage of Data Known")
ylabel("Relative Frobenious Norm Error")
legend(["Trace Inverse","Nuclear Norm", "Noisy matrix"])
set(gca,'fontsize', 12);

error
