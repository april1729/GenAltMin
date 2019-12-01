clear
n=100;
m=200;
r=10;
noiseLevel=.1;
numRuns=10;


p=0.4;
mu_range=[0.1,1,10,100,1000];
[D, A, b,M]=generateMatrixCompletionProblem(m,n,r,p, noiseLevel);

for run=1:1
    figure()
    for i=1:length(mu_range)
        
        
        
        opts.r=20;
        opts.maxIter=200;
        opts.mu=mu_range(i);
        opts.xTol=1e-5;
        opts.f=@(x) 100*50./(50+x).^2;
        
        opts.obj=@(U,V) obj(U*V', @(x) 5*x, opts.f, opts.mu, A, b);
        
        [ U_ti,V_ti , obj_nuc] = GenAltMinSD(M,A,b,opts );
        error(i, 1,run)=norm(U_ti*V_ti'-D,'fro')/norm(D, 'fro');
        rank_list(i,1,run)=sum(svd(U_ti*V_ti')>0.0001);
        epsilon(i,1,run)=norm(A*vec(U_ti*V_ti')-b,'fro');
        
        opts.f=@(x) 1;
        
        [ U_nuc, V_nuc, obj_ti] = GenAltMinSD(M,A,b,opts );
        error(i, 2, run)=norm(U_nuc*V_nuc'-D,'fro')/norm(D, 'fro');
        rank_list(i,2,run)=sum(svd(U_nuc*V_nuc')>0.0001);
        epsilon(i,2,run)=norm(A*vec(U_nuc*V_nuc')-b,'fro');
        subplot(length(mu_range),2,2*i-1)
        bar([svds(D,10), svds(U_nuc*V_nuc',10), svds(U_ti*V_ti',10)])
        xlabel("Index")
        ylabel("Singular Value")
        title("\mu="+opts.mu)
        %legend(["Original Matrix", "Nuclear Norm", "Trace Inverse"])
        set(gca,'fontsize', 10);
        
        subplot(length(mu_range),2,2*i)
        data1=svds(D,20);
        data2=svds(U_nuc*V_nuc',20);
        data3=svds(U_ti*V_ti',20);
        bar(11:20,[data1(11:20),data2(11:20),data3(11:20) ])
        %xlabel("Index")
        title("\mu="+opts.mu)
        % ylabel("Singular Value")
        
        legend(["Original Matrix", "Nuclear Norm", "Trace Inverse"])
        set(gca,'fontsize', 10);
        
    end
end

figure()
semilogx(epsilon, mean(error,3),'linewidth', 4)
hold on
eps_range=min(min(epsilon)):0.1:max(max(epsilon));
semilogx(eps_range, ones(size(eps_range))*norm(A*vec(D)-b)/norm(b),'--','linewidth', 4)
xlabel("\epsilon")
ylabel("Relative Frobenious Norm Error")
legend(["Trace Inverse","Nuclear Norm", "Noisy matrix"])
set(gca,'fontsize', 24);