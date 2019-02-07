%% Sparse Optimization MM Algorithm
% I coded a funtion sparseOptimizationMM.m that does an MM algorithm for
% for any regularizer, so far i have SCAD, LogDet, and the LPCC
% regularizers programmed for it.  The code can acomidate any regularizer
% with just a couple extra lines of code, and takes the name of the
% regularizer as input.
%
% Its good code.

m=200;
n=10000;
s=30;

A=(rand(m,n)-0.5)*10;
xSoln=zeros(n,1);
xSoln(datasample(1:n, s))=0.5+rand(s,1)*5;
b=A*xSoln;

[xSCAD,obj]=sparseOptimizationMM(zeros(n,1), A,b, {"SCAD", 0.5,2});
[xLogDet,obj]=sparseOptimizationMM(zeros(n,1), A,b, {"LogDet",0.5 });
[x,obj]=sparseOptimizationMM(zeros(n,1), A,b, {"MR", 0.1});
xL1 = linprog(ones(n,1),[],[],A,b, zeros(n,1), []);

fprintf("SCAD sparsity:   %i\n", sum(xSCAD>0));
fprintf("LogDet sparsity: %i\n", sum(xSCAD>0));
fprintf("LPCC sparsity:   %i\n", sum(x>0));
fprintf("L1 sparsity:     %i\n", sum(xL1>0));



%% A comparison of all these methods

n=10000;
s=30;
spasity=[];
error=[];
mList=100:10:250;
if exist("sparse_mm_comparison.mat", 'file')~=2
for i=1:length(mList)
    m=mList(i);
    A=(rand(m,n)-0.5)*10;
    xSoln=zeros(n,1);
    xSoln(datasample(1:n, s))=rand(s,1)*5;
    b=A*xSoln;
    
    [xSCAD,obj]=sparseOptimizationMM(zeros(n,1), A,b, {"SCAD", 0.075,2.5});
    [xLogDet,obj]=sparseOptimizationMM(zeros(n,1), A,b, {"LogDet",0.1 });
    [xLPCC,obj]=sparseOptimizationMM(zeros(n,1), A,b, {"MR", 0.1});
    xL1 = linprog(ones(n,1),[],[],A,b, zeros(n,1), []);
    
    
    sparity(1,i)= sum(xSCAD~=0);
    sparity(2,i)= sum(xLogDet~=0);
    sparity(3,i)= sum(xLPCC~=0);
    sparity(4,i)= sum(xL1~=0);
    
    error(1,i)=norm(xSCAD-xSoln)/norm(xSoln);
    error(2,i)=norm(xLogDet-xSoln)/norm(xSoln);
    error(3,i)=norm(xLPCC-xSoln)/norm(xSoln);
    error(4,i)=norm(xL1-xSoln)/norm(xSoln);
    
end
save("sparse_mm_comparison.mat")
else
    load("sparse_mm_comparison.mat")
end

[rankTable, errorTable]=displayResults(mList, sparity', error', [], {'Number_of_Constraints', 'SCAD', 'LogDet', 'LPCC', 'l1_Norm'});

%% Extending to PSD Rank Minimization
% I now use pretty much the same algorithm for rank mimization.  In
% general, for a regularizer $f(x_i)$, the weight matrix $W^n$ is defined
% as
% $$W^n={P^n}^T diag(f'(\sigma_1), f'(\sigma_2), ..., f'(\sigma_n)) P^n$$
% Where
% $$X^n={P^n}^T diag(\sigma_1, \sigma_2, ..., \sigma_n) P^n$$

n=50;
r=3;
p=0.2;
rankTol=0.00001;
[D, A, b]=generateMatrixCompletionProblemPSD(n,r,p,0);
[XMR,obj]=rankMinimizationPSDMM(zeros(n,n), [], [], 0, A, b, {"MR",1});
[XlogDet,obj]=rankMinimizationPSDMM(zeros(n,n), [], [], 0, A, b, {"LogDet",1});
[XSCAD,obj]=rankMinimizationPSDMM(zeros(n,n), [], [], 0, A, b, {"SCAD",4,3});
[ Xnuc ] = nuclearNormPSD(zeros(n,n),[],[],0,A,b, []);

fprintf("SCAD rank:   %i\n", sum(svd(XSCAD)>rankTol));
fprintf("MR rank:     %i\n", sum(svd(XMR)>rankTol));
fprintf("LogDet rank: %i\n", sum(svd(XlogDet)>rankTol));
fprintf("NNM rank:    %i\n", sum(svd(Xnuc)>rankTol));


fprintf("SCAD error:   %f\n", norm(XSCAD-D,'fro')/norm(D,'fro'));
fprintf("MR error:     %f\n", norm(XMR-D,'fro')/norm(D,'fro'));
fprintf("LogDet error: %f\n", norm(XlogDet-D,'fro')/norm(D,'fro'));
fprintf("NNM error:    %f\n", norm(Xnuc-D,'fro')/norm(D,'fro'));


figure()
hold on
plot(svd(Xnuc), '-')
plot(svd(XlogDet),'--')
plot(svd(XMR),':')
plot(svd(XSCAD), '-.')
%% Nonsymmetric Rank Minimization
% Again, not much really changes here.  Replace $X$ with $X^TX$ everywhere
% besides the constraints.

m=100;
n=50;
r=3;
p=0.2;
rankTol=0.1;
X0=zeros(m,n);
[D, A, b]=generateMatrixCompletionProblem(n,m,r,p,0);

[ Xnuc ] = nuclearNorm(X0,[],[],0,A,b, []);
[XMR,obj]=rankMinimizationMM(Xnuc, [], [], 0, A, b, {"MR",3});
[XlogDet,obj]=rankMinimizationMM(Xnuc, [], [], 0, A, b, {"LogDet",3});
[XSCAD,obj]=rankMinimizationMM(Xnuc, [], [], 0, A, b, {"SCAD",2,4});


fprintf("SCAD rank:   %i\n", sum(svd(XSCAD)>rankTol));
fprintf("MR rank:     %i\n", sum(svd(XMR)>rankTol));
fprintf("LogDet rank: %i\n", sum(svd(XlogDet)>rankTol));
fprintf("NNM rank:    %i\n", sum(svd(Xnuc)>rankTol));

fprintf("SCAD error:   %f\n", norm(XSCAD-D,'fro')/norm(D,'fro'));
fprintf("MR error:     %f\n", norm(XMR-D,'fro')/norm(D,'fro'));
fprintf("LogDet error: %f\n", norm(XlogDet-D,'fro')/norm(D,'fro'));
fprintf("NNM error:    %f\n", norm(Xnuc-D,'fro')/norm(D,'fro'));


figure()
hold on
plot(svd(Xnuc), '-')
plot(svd(XlogDet),'--')
plot(svd(XMR),':')
plot(svd(XSCAD), '-.')

%% Sparsity Regulated Rank Minimzation
%
% now lets do an MM algorithm that combines both of these.  For a surrogate
% function, we just need to add the surrogate functions for rank and l0.
% we can mix and match the different regularizers if we want, i suppose.

[D,S, A, b]=generateRobustMatrixCompletionProblem(100,50,5,500,0.7, 0);
[ Xc, Sc ] = sparsityRegulatedRankMinimizationConvex( zeros(100,50), [],[],0, A,b, {0.1});
[ Xmm, Smm ] = sparsityRegulatedRankMinimizationMM( Xc, [],[],0, A,b, {0.1,{'lpcc', 0.5},{'LogDet', 2}});
sum(sum(abs(full(Sc))>0.00001))
sum(sum(abs(full(Smm))>0.00001))
norm(Xc-D, 'fro')/norm(D,'fro')
norm(Xmm-D, 'fro')/norm(D,'fro')

figure()
heatmap(abs(Sc))
figure
heatmap(abs(Smm))