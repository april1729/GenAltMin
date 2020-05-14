load_path

[M_train_all, M_test_all] = loadJester();
rank = 10;

for fold=1:5
	opts=[];
	M_test=M_test_all{fold};
	M_train=M_train_all{fold};

	n_samples=sum(sum(M_test~=0));
	[n_movies,n_users]=size(M_train);

	[row_test, col_test]=find(M_test);
	b_test=vec(M_test);
	b_test=b_test(sub2ind(size(M_test), row_test, col_test));
	SDopts.maxIter=1;
	opts.SDopts=SDopts;

	opts.r=rank;
	opts.obj=@(U,V) sum(abs(vec(sparse_multiply(U, V, row_test, col_test, n_movies,n_users))-vec(M_test)))/n_samples;
	opts.beta=1e-6;
	opts.fTol=0;
	opts.maxIter=40;
	opts.gamma=10;
	[ Ut,Vt , obj_ti] = GenASD(M_train,opts );
	M_predict_t1=sparse_multiply(Ut, Vt, row_test, col_test, n_movies,n_users);

	b_predict_ti=vec(M_predict_t1);
	b_predict_ti=b_predict_ti(sub2ind(size(M_test), row_test, col_test));
	b_predict_ti=(min(max(b_predict_ti, -10),10));



	opts.r=rank;
	opts.obj=@(U,V) sum(abs(vec(sparse_multiply(U, V, row_test, col_test, n_movies,n_users))-vec(M_test)))/n_samples;
	opts.beta=1;
	opts.fTol=0;
	opts.maxIter=40;
	opts.gamma=10;
	opts.f=@(x,gamma) 1;
	[ U,V , obj_nn] = GenASD(M_train,opts );
	M_predict=sparse_multiply(U, V, row_test, col_test, n_movies,n_users);
	svd_ti=svd(U'*U+V'*V);
	b_predict_nn=vec(M_predict);
	b_predict_nn=b_predict_nn(sub2ind(size(M_test), row_test, col_test));
	b_predict_nn=(min(max(b_predict_nn,-10),10));




	optsLMA.est_rank = 1;
	optsLMA.print=2;
	Idx=find(M_train~=0);
	MIdx= M_train(Idx);
	tic
	[X,Y] = lmafit_mc_adp(size(M_train,1),size(M_train,2),rank,Idx,MIdx,optsLMA);
	time=toc;
	M_predict_LMa=sparse_multiply(X, Y', row_test, col_test, n_movies,n_users);

	b_predict_LMaFit=vec(M_predict_LMa);
	b_predict_LMaFit=b_predict_LMaFit(sub2ind(size(M_test), row_test, col_test));
	b_predict_LMaFit=(min(max(b_predict_LMaFit, -10),10));
	nmae(fold,1)=mean(abs(full(b_test-b_predict_ti)))/20;
	nmae(fold,2)=mean(abs(full(b_test-b_predict_nn)))/20;
	nmae(fold,3)=mean(abs(full(b_test-b_predict_LMaFit)))/20
	%fprintf("Trace Inverse MAE: %f \n",mean(abs(full(b_test-b_predict_ti)))/20);
	%fprintf("Nucler Norm MAE:   %f \n",mean(abs(full(b_test-b_predict_nn)))/20);
	%fprintf("LMaFit MAE:        %f \n",mean(abs(full(b_test-b_predict_LMaFit)))/20);
end

    nmae(6,1)=mean(nmae(1:5,1));
    nmae(6,2)=mean(nmae(1:5,2));
    nmae(6,3)=mean(nmae(1:5,3));
    nmae=round(nmae, 4)
