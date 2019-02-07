if exist("nov28Section1.mat", 'file')~=2    
    plist=0.3:0.05:0.6;
    rank=[];
    error=[];
    
    rankTol=.1;

    for i=1:length(plist)
        [D, A, b]=generateMatrixCompletionProblem(15,3,plist(i),0);
        
        [ Xnuc ] = nuclearNormPSD(zeros(15,15),[],[],0,A,b, []);
        [Xar, U, obj]= aggregateRelaxation(zeros(15,15),[],[],0,A,b,[10,0.9] );

        error(i,1)=norm(Xnuc-D, 'fro')/norm(D, 'fro');
        error(i,2)=norm(Xar-D, 'fro')/norm(D, 'fro');
        
        rank(i,1)=sum(svd(Xnuc)>rankTol);
        rank(i,2)=sum(svd(Xar)>rankTol);
                
    end
    
    save("nov28Section1.mat")
else
    load("nov28Section1.mat")
end
[rankTable, errorTable, runTimeTable]=displayResults(plist, rank, error, [], {'Percent_of_Data_Missing', 'Nuclear_norm', 'aggregate_relaxation'});