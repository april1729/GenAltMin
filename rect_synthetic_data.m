clear
n=100;
m=200;
r=10;
noiseLevel=.1;
numRuns=10;


p_range=0.2:0.1:1;
for run=1:1
    figure()
    for i=1:length(p_range)
        [D, A, b,M]=generateMatrixCompletionProblem(m,n,r,p_range(i), noiseLevel);
        
        
        
        opts.r=20;
        opts.maxIter=2000;
        opts.mu=1;
        opts.xTol=1e-7;
        opts.f=@(x) 100*50./(50+x).^2;
        
        opts.obj=@(U,V) obj(U*V', @(x) 5*x, opts.f, opts.mu, A, b);
        
        [ U_ti,V_ti , obj_nuc] = GenASD(M,A,b,opts );
        error(i, 1,run)=norm(U_ti*V_ti'-D,'fro')/norm(D, 'fro');
        rank_list(i,1,run)=sum(svd(U_ti*V_ti')>0.0001);
        
        
        opts.f=@(x) 1;
        
        [ U_nuc, V_nuc, obj_ti] = GenASD(M,A,b,opts );
        error(i, 2, run)=norm(U_nuc*V_nuc'-D,'fro')/norm(D, 'fro');
        rank_list(i,2,run)=sum(svd(U_nuc*V_nuc')>0.0001);
        
        subplot(length(p_range),2,2*i-1)
        bar([svds(D,10), svds(U_nuc*V_nuc',10), svds(U_ti*V_ti',10)])
        xlabel("Index")
        ylabel("Singular Value")
        title(p_range(i)+"%")
        %legend(["Original Matrix", "Nuclear Norm", "Trace Inverse"])
        set(gca,'fontsize', 10);
        
        subplot(length(p_range),2,2*i)
        data1=svds(D,20);
        data2=svds(U_nuc*V_nuc',20);
        data3=svds(U_ti*V_ti',20);
        bar(11:20,[data1(11:20),data2(11:20),data3(11:20) ])
        %xlabel("Index")
        title(p_range(i)+"%")
       % ylabel("Singular Value")
        
        %legend(["Original Matrix", "Nuclear Norm", "Trace Inverse"])
        set(gca,'fontsize', 10);
        
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


