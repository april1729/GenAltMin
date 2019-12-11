clear
load test_parameters
T=table('RowNames', {'m', 'n', 'r', 'noise', 'p',  '__','GenAltMin (Trace Inverse)', 'GenASD (Trace Inverse)', 'FaNCL',...
     'SVT', 'LMaFit (correct r)', 'LMaFit (incorrect r)'});
temp=zeros(5, 2*length(test_parameters));
temp(1:5, 1:2:2*length(test_parameters)-1)=(test_parameters(1:5,:));
temp(1:5, 2:2:2*length(test_parameters))=(test_parameters(1:5,:));
T(1:5, 1:2*length(test_parameters))=num2cell(temp);
%T(6,1:2*length(test_parameters))=cellstr(repmat(["RFNE","Time"], 1,8));


%for i=1:length(test_parameters)
    i=1;
    m=300;%test_parameters(1,i);
    n=200;%test_parameters(2,i);
    r=5;%test_parameters(3,i);
    noise=0.1;%test_parameters(4,i);
    p=0.4;%test_parameters(5,i);
    
    r_upper_bound=2*r;
    beta_genASD=0.01;%test_parameters(6,i);
    beta_genAltMin=5e-4;%test_paraclearmeters(7,i);
    beta_nuclear_norm=1;%test_parameters(8,i);
    
    
    
    [D, A, b,M]=generateMatrixCompletionProblem(m,n,r,p, noise);
    [~,j] = find(A);
    [rows,cols]=ind2sub([m,n], j);
    Omega=(M~=0);
    opts.obj=@(U,V) sum(1- gamma./(gamma+eig([U;V]'*[U;V])))...
        +(opts.mu/2)*norm(A*vec(sparse_multiply(U,V, rows, cols,m,n))-b)^2;
    
    
    % GenAltMin (Trace Inverse)
    gamma=norm(b)/sqrt(2*p*r_upper_bound);
    
    opts.r=r_upper_bound;
    opts.f=@(x) gamma./(gamma+x).^2;
    opts.mu=beta_genAltMin;
    opts.obj=@(U,V) sum(1- gamma./(gamma+eig([U;V]'*[U;V])))...
        +(opts.mu/2)*norm(A*vec(sparse_multiply(U,V, rows, cols,m,n))-b)^2;
    
    tic;
    [ U,V , obj,~] = genAltMin_v2(M,Omega,opts );
    time=toc;
    RFNE=norm(D-U*V', 'fro')/norm(D, 'fro');
    T(7,2*i-1)={RFNE};
    T(7,2*i)={time};
    % GenASD (Trace Inverse)
    opts.mu=beta_genASD;
    
    tic;
    [ U,V , obj,~] = GenASD(M,A,b,opts );
    time=toc;
    RFNE=norm(D-U*V', 'fro')/norm(D, 'fro');
    
    T(8,2*i-1)={RFNE};
    T(8,2*i)={time}
    
    % gSVT
    
    
    tic;
    para.maxR=r_upper_bound;
    para.maxIter=1000;
    para.tol=1e-4;
    para.regType=1;
    [U, S, V, output ] = FastMatComp( M, 1/beta_genASD, gamma, para);
    time=toc;
    X=U*S*V';
    RFNE=norm(D-X,'fro')/norm(D,'fro');
    
    T(9,2*i-1)={RFNE};
    T(9,2*i)={time}
    
    
    % GenAltMin (Nuclear Norm)
    
%     opts.f=@(x) 1;
%     opts.mu=beta_nuclear_norm;
%     
%     opts.obj=@(U,V) sum(eig([U;V]'*[U;V]))...
%         +(opts.mu/2)*norm(A*vec(sparse_multiply(U,V, rows, cols,m,n))-b)^2;
%     
%     tic;
%     [ U,V , ~,~] = GenAltMin(M,A,b,opts );
%     time=toc;
%     RFNE=norm(D-U*V', 'fro')/norm(D, 'fro');
%     T(10,2*i-1)={RFNE};
%     T(10,2*i)={time}
%     
%     % GenASD (Nuclear Norm)
%     
%     
%     tic;
%     [ U,V , ~,~] = GenASD(M,A,b,opts );
%     time=toc;
%     RFNE=norm(D-U*V', 'fro')/norm(D, 'fro');
%     T(11,2*i-1)={RFNE};
%     T(11,2*i)={time}
%     
    % SVT
    tic;
    [U,S,V,numiter]  = FPC([m,n],find(M~=0),b,1/beta_nuclear_norm);
    time=toc;
    RFNE=norm(D-U*S*V','fro')/norm(D,'fro');
    T(12,2*i-1)={RFNE};
    T(12,2*i)={time}
    
    
    % LMaFit (correct r)
    optsLMA.est_rank = 1;
    Idx=find(M~=0);
    MIdx= M(Idx);
    rank = 2*r;
    tic
    [X,Y] = lmafit_mc_adp(size(M,1),size(M,2),rank,Idx,MIdx,optsLMA);
    time=toc;
    RFNE=norm(D-X*Y, 'fro')/norm(D, 'fro');
    T(13,2*i-1)={RFNE};
    T(13,2*i)={time}
    
    % LMaFit (incorrect r)
    rank = 1.5*r;
    tic
    [X,Y] = lmafit_mc_adp(size(M,1),size(M,2),rank,Idx,MIdx,optsLMA);
    time=toc;
    RFNE=norm(D-X*Y, 'fro')/norm(D, 'fro');
    T(14,2*i-1)={RFNE};
    T(14,2*i)={time}
    
    
%end