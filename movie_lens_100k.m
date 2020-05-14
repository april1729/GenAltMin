load_path

for fold=1:5
    [A,b,M_train, M_test] = loadMovieLens100k(fold);
    [n_movies,n_users]=size(M_train)
    [row_test, col_test]=find(M_test);
    b_test=vec(M_test);
    b_test=b_test(sub2ind(size(M_test), row_test, col_test));
    SDopts.maxIter=1;
    opts.SDopts=SDopts;
    
    opts.r=10;
    opts.obj=@(U,V) sum(abs(vec(sparse_multiply(U, V, row_test, col_test, 943, 1682))-vec(M_test)))/20000;
    opts.beta=5e-3;
    opts.fTol=1e-6;
    opts.maxIter=100;
    opts.gamma=10;
    [ U,V , obj_ti] = GenASD(M_train,opts );
    M_predict_t1=sparse_multiply(U, V, row_test, col_test, 943, 1682);
    svd_ti=svd(U'*U+V'*V);
    b_predict_ti=vec(M_predict_t1);
    b_predict_ti=b_predict_ti(sub2ind(size(M_test), row_test, col_test));
    b_predict_ti=round(min(max(b_predict_ti, 1),5));
    
    [U,S,V,numiter]  = FPC([n_movies,n_users],find(M_train~=0),M_train(find(M_train~=0)),1);

        
    % opts.r=10;
    % opts.obj=@(U,V) sum(abs(vec(sparse_multiply(U, V, row_test, col_test, 943, 1682))-vec(M_test)))/20000;
    % opts.beta=1;
    % opts.fTol=1e-6;
    % opts.maxIter=100;
    % opts.gamma=10;
    % opts.f=@(x,gamma) 1;
    % [ U,V , obj_ti] = GenASD(M_train,opts );
     M_predict=sparse_multiply(U, V, row_test, col_test, 943, 1682);
    % svd_ti=svd(U'*U+V'*V);
     b_predict_nn=vec(M_predict);
     b_predict_nn=b_predict_nn(sub2ind(size(M_test), row_test, col_test));
     b_predict_nn=round(min(max(b_predict_nn, 1),5));
     opts=[];

    
    
    optsLMA.est_rank = 1;
    optsLMA.print=0;
    Idx=find(M_train~=0);
    MIdx= M_train(Idx);
    rank = 1;
    tic
    [X,Y] = lmafit_mc_adp(size(M_train,1),size(M_train,2),rank,Idx,MIdx,optsLMA);
    time=toc;
    M_predict_LMa=sparse_multiply(X, Y', row_test, col_test, 943, 1682);
    
    b_predict_LMaFit=vec(M_predict_LMa);
    b_predict_LMaFit=b_predict_LMaFit(sub2ind(size(M_test), row_test, col_test));
    b_predict_LMaFit=round(min(max(b_predict_LMaFit, 1),5));
    
    MAE_Table(fold,1)=fold;
    MAE_Table(fold,2)=mean(abs(full(b_test-b_predict_ti)))/4;
    MAE_Table(fold,3)=mean(abs(full(b_test-b_predict_nn)))/4;

    MAE_Table(fold,4)=mean(abs(full(b_test-b_predict_LMaFit)))/4
    
    
end
    MAE_Table(6,2)=mean(MAE_Table(1:5,2));
    MAE_Table(6,3)=mean(MAE_Table(1:5,3));
    MAE_Table(6,4)=mean(MAE_Table(1:5,4));
    MAE_Table=round(MAE_Table, 4)
