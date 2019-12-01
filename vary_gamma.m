n=300;
m=200;
r=5;
p=0.3;
noiseLevel=.1;
numRuns=10;
[D, A, b,M]=generateMatrixCompletionProblem(m,n,r,p, noiseLevel);


gamma_range=[0.03125, 0.0625,0.125, 0.25, 0.5,1,2,4,8,16,32,64,128, 256];
for run=1:1
    figure()
    opts.r=10;
    opts.maxIter=200;
    opts.mu=0.01;
    opts.xTol=1e-5;
    opts.f=@(x) 1/norm(D);
    opts.obj=@(u,v) 0;
    [ U_n,V_n , obj_nuc] = GenASD(M,A,b,opts );
    
    for i=1:length(gamma_range)
        
        
        
        opts.r=20;
        opts.maxIter=200;
        opts.mu=0.01;        

        opts.xTol=1e-5;
        opts.f=@(x) gamma_range(i)./(gamma_range(i)+x).^2;
        opts.obj=@(U,V) obj(U*V', @(x) 5*x, opts.f, opts.mu, A, b);
        
        [ U_ti,V_ti , obj_nuc] = GenASD(M,A,b,opts );
        error(i, 1,run)=norm(U_ti*V_ti'-D,'fro')/norm(D, 'fro')
        rank_list(i,1,run)=sum(svd(U_ti*V_ti')>0.0001);
    end
end

figure()
semilogx(gamma_range, error,'-o','linewidth',  2)
hold on
semilogx(gamma_range, ones(size(gamma_range))*norm(U_n*V_n'-D,'fro')/norm(D, 'fro'),'--','linewidth', 2)
semilogx(gamma_range, ones(size(gamma_range))*norm(A*vec(D)-b)/norm(b),'--','linewidth', 2)

xlabel("\gamma")
ylabel("Relative Frobenious Norm Error")
legend(["Trace Inverse","Nuclear Norm", "Noisy matrix"])
set(gca,'fontsize', 12);