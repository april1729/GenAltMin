n=1000;        opts.mu=mu_range(i);

m=500;
r=5;
p=0.1;
noiseLevel=.05;
numRuns=10;
[D, A, b,M]=generateMatrixCompletionProblem(m,n,r,p, noiseLevel);


gamma_range=[2,4,8,16,32,64,128, 256,512,1024,2048,2^12,2^13 ];
for run=1:1
    figure()
    opts.r=10;
    opts.maxIter=200;
    opts.mu=0.01;
    opts.xTol=1e-5;
    opts.f=@(x, gamma) 1/norm(D);
    opts.obj=@(u,v) 0;
    [ U_n,V_n , obj_nuc] = GenASD(M,opts );
    
    for i=1:length(gamma_range)
        
        
        
        opts.r=10;
        opts.maxIter=1000;
        opts.mu=0.01;        

        opts.xTol=1e-5;
        opts.f=@(x, gamma) gamma_range(i)./(gamma_range(i)+x).^2;
        opts.obj=@(U,V) obj(U*V', @(x) 5*x, opts.f, opts.mu, A, b);
        opts.gamma=gamma_range(i);
        opts.gamma0=100*gamma_range(i);
        opts.beta=0.001;
        opts.beta0=0.001;
        [ U_ti,V_ti , obj_nuc] = genAltMin_v2(M,(M~=0),opts );
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
ylabel("Relative Frobenious Norm Error")
legend(["Trace Inverse","Nuclear Norm", "Noise Matrix"])
set(gca,'fontsize', 12);