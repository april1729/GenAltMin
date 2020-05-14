load_path
n=1000;        
m=500;
r=5;
p=0.1;
noiseLevel=.05;

numRuns=10;



[D, A, b,M]=generateMatrixCompletionProblem(m,n,r,p, noiseLevel);


gamma_range=2.^(-2:1:8);
for run=1:numRuns
    figure()
    opts.r=10;
    opts.maxIter=200;
    opts.mu=0.000001;
    opts.xTol=1e-9;
    opts.f=@(x, gamma) 1/norm(D);
    opts.obj= @(U,V) obj(U*V', @(x)1/norm(D) , opts.f, opts.mu, A, b);;
    [ U_n,V_n , obj_nuc] = GenASD(M,opts );
    
    for i=1:length(gamma_range)
        
        
        
        opts.r=10;
        opts.maxIter=200;
        opts.mu=0.01;        

        opts.xTol=1e-5;
        opts.f=@(x, gamma) gamma_range(i)./(gamma_range(i)+x).^2;
        opts.obj=@(U,V) obj(U*V', @(x) 5*x, opts.f, opts.mu, A, b);
        opts.gamma=gamma_range(i);
        opts.gamma0=gamma_range(i);
        opts.beta=0.01;
        opts.beta0=0.01;
        [ U_ti,V_ti , obj_nuc] = GenASD(M,opts );
        error(i, 1,run)=norm(U_ti*V_ti'-D,'fro')/norm(D, 'fro')
        rank_list(i,1,run)=sum(svd(U_ti*V_ti')>0.0001);
    end
end

figure()
semilogx(gamma_range, mean(error,3),'-o','linewidth',  2)
hold on
semilogx(gamma_range, ones(size(gamma_range))*norm(U_n*V_n'-D,'fro')/norm(D, 'fro'),'--','linewidth', 2)
semilogx(gamma_range, ones(size(gamma_range))*norm(A*vec(D)-b)/norm(b),'--','linewidth', 2)

xlabel("\gamma")
ylabel("RFNE")
legend(["Trace Inverse","Nuclear Norm", "Noise Matrix"])
set(gca,'fontsize', 12);