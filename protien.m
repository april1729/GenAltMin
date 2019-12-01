opts.exact=false;
opts.maxIter=70;
opts.mu=1e-4;
opts.r=10;
opts.xtol=5e-5;
files=[ "1AX8.mat", "1KDH.mat", "1HOE.mat", "1BPM.mat","1PTQ.mat", "1RGS.mat"]; %",
percents=[100];
output_table=[];
for file=files
    for percent=percents
        load(file)
        %opts.obj=@(x) norm(calcDistance(x*x')-prob.Dmat100,'fro')/norm(prob.Dmat100,'fro');
        
        P0=rand(length(prob.A),10);
        %[A,b]=constructProteinProblem(file, 100);
        D=calcDistance(prob.A'*prob.A);
        noise_matrix=randn(size(D));
        noise_matrix=noise_matrix*noise_matrix';
        [A, b] = sampleUniformGram( D+0.2*noise_matrix,0.05);
        new_row(1)=file;
        opts.mu=1e-4;
        opts.f=@(x) 100/(100+x).^2;
        tic;
        [X,obj]=gnrtn(P0,A,b,opts);
        
        opts.mu=1e-3;
        new_row(3)=toc;
        
        new_row(2)=norm(calcDistance(X*X')-D,'fro')/norm(D,'fro');
        
        opts.f=@(x) 1;
        
        tic;
        [X,obj]=gnrtn(P0,A,b,opts);
        new_row(5)=toc;
        new_row(4)=norm(calcDistance(X*X')-D,'fro')/norm(D,'fro')
        
        output_table=[output_table;new_row];
    end
end
