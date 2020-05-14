load_path
clear
m=1000;
n=500;
r=10;
noiseLevel=.1;
numRuns=20;

p=0.1;


opts.r=2*r;
opts.maxIter=2000;
opts.xTol=1e-7;
opts.beta=0.1;
gamma=50;
opts.gamma=gamma;



beta_range=[1,0.1,0.05,0.01,0.001,1e-4]
for run=1:numRuns
    figure()
    
    [D, A, b,M]=generateMatrixCompletionProblem(m,n,r,p, noiseLevel);
    Omega=(M~=0);
    [~,j] = find(A);
    [rows,cols]=ind2sub([m,n], j);
    opts.obj=@(U,V) sum(1- gamma./(gamma+eig([U;V]'*[U;V])))...
        +(0.001/2)*norm(A*vec(sparse_multiply(U,V, rows, cols,m,n))-b)^2;
    
    for i=1:length(beta_range)
        beta=beta_range(i);
        opts.beta=beta;
        
        
        [ U_ti,V_ti , obj_nuc] = GenASD(M,opts );
        error(i, 1,run)=norm(U_ti*V_ti'-D,'fro')/norm(D, 'fro');
        rank_list(i,1,run)=sum(svd(U_ti*V_ti')>0.0001);
        
        [X]=svt(M, Omega,opts.beta, 2*r);
        
        error(i, 2, run)=norm(X-D,'fro')/norm(D, 'fro');
        rank_list(i,2,run)=sum(svd(X)>0.0001);
        
        subplot(length(beta_range),2,2*i-1)
        bar([svds(D,r), svds(X,r), svds(U_ti*V_ti',r)])
        xlabel("Index")
        ylabel("Singular Value")
        title("\beta="+beta)
        %legend(["Original Matrix", "Nuclear Norm", "Trace Inverse"])
        set(gca,'fontsize', 10);
        
        subplot(length(beta_range),2,2*i)
        data1=svds(D,20);
        data2=svds(X,2*r);
        data3=svds(U_ti*V_ti',2*r);
        bar(r+1:2*r,[data1(r+1:2*r),data2(r+1:2*r),data3(r+1:2*r) ])
        %xlabel("Index")
        % ylabel("Singular Value")
        
        %
        set(gca,'fontsize', 10);
        
    end
end
legend(["Original Matrix", "Nuclear Norm", "Trace Inverse"])

figure()
semilogx(beta_range, mean(error,3),'linewidth', 4)
xlabel("\beta")
ylabel("Relative Frobenious Norm Error")
legend(["Trace Inverse","Nuclear Norm", "Noisy matrix"])
set(gca,'fontsize', 12);

error
