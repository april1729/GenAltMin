%% December 4th
%% Testing nuclear norm constraint

A=vec(eye(10));
b=1;
[X,obj]=matrixRelaxationMM(zeros(10,10), [], [], 0, A', b, [1,1]);
sparse(X)



[X,obj]=matrixRelaxationMM(zeros(10,10), [], [], 0, A', 10, [0.001,1]);
sparse(X)
%%
% Probably just a numerical percision thing?

%% KKT points of this problem:


%% Testing aggregate relaxation more:
if exist("dec4Section1.mat", 'file')~=2
    rankTol=0.01;
    plist=0.1:0.1:0.5;
    for i =1:length(plist)
        [D, A, b]=generateMatrixCompletionProblem(50,3,plist(i),0);
        
        [ Xnuc ] = nuclearNormPSD(zeros(50,50),[],[],0,A,b, []);
        [Xar, U, obj]= aggregateRelaxation(zeros(50,50),[],[],0,A,b,[20,0.9] );
        [Xmm, obj]= matrixRelaxationMM(zeros(50,50),[],[],0,A,b,[5,1] );
        
        error(i,1)=norm(Xnuc-D, 'fro')/norm(D, 'fro');
        error(i,2)=norm(Xar-D, 'fro')/norm(D, 'fro');
        error(i,3)=norm(Xmm-D, 'fro')/norm(D, 'fro');
        
        rank(i,1)=sum(svd(Xnuc)>rankTol);
        rank(i,2)=sum(svd(Xar)>rankTol);
        rank(i,3)=sum(svd(Xmm)>rankTol);
        
    end
    
    save("dec4Section1.mat")
else
    load("dec4Section1.mat")
end
[rankTable, errorTable, runTimeTable]=displayResults(plist, rank, error, [], {'Percent_of_Data_Missing', 'Nuclear_norm', 'Aggregate_Relaxation', 'Matrix_Relaxation'});