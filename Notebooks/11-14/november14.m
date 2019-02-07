%% Testing SeDuMi for Aggregated Relaxation Algorithm
%
% Lets set up a test case real quick to test this on.
clear;
n=100;
r=30;

[points, distances, A, b]=generateSensorProblem(n,r,.5, 0);
U=eye(n);
eps=1;
c=10;
mu=1;

tic;
cvx_begin
cvx_solver sedumi
variable X(n,n)
minimize mu*pos(sum(sum(X.*U))-eps)+(c/2)* square_pos(sum(sum(X.*U))-eps)
subject to
A*vec(X)==b;
X == semidefinite(n)
cvx_end
cvx_time=toc;


tic;
cvx_begin
cvx_solver sdpt3
variable Xsdpt3cvx(n,n)
minimize mu*pos(sum(sum(Xsdpt3cvx.*U))-eps)+(c/2)* square_pos(sum(sum(Xsdpt3cvx.*U))-eps)
subject to
A*vec(Xsdpt3cvx)==b;
Xsdpt3cvx == semidefinite(n)
cvx_end
sdpt3cvx_time=toc;
%%
% Now we solve this with sedumi
numCon=length(b);

Asedumi=zeros(numCon+3, n^2+5);
Asedumi(1:numCon,6:n^2+5)=A;
Asedumi(numCon+1,6:n^2+5)=vec(U)';
Asedumi(numCon+1:numCon+3, 1:5)=[-1, 1, 0, 0, 0;1, 0, 0 ,0 -1;0,0,0,1,0];
csedumi=zeros(5+n^2, 1);
csedumi(1)=mu;
csedumi(3)=c/2;
bsedumi=[b;eps;0;0.5];

K.l=2;
K.r=3;
K.s=n;
tic;
[x,y,info] = sedumi(Asedumi',bsedumi,csedumi,K);
sedumi_time=toc;
Xsedumi=mat(x(6:n^2+5));

%%
% Well, it worked, lets look at the results then.


fprintf("Sedumi Time:          %f \t Objective: %f \n", sedumi_time, mu*pos(sum(sum(Xsedumi.*U))-eps)+(c/2)* pos(sum(sum(Xsedumi.*U))-eps)^2);
fprintf("CVX with Sedumi Time: %f \t Objective: %f \n", cvx_time, mu*pos(sum(sum(X.*U))-eps)+(c/2)* pos(sum(sum(X.*U))-eps)^2);
fprintf("CVX with SDPT3  Time: %f \t Objective: %f \n", sdpt3cvx_time, mu*pos(sum(sum(Xsdpt3cvx.*U))-eps)+(c/2)* pos(sum(sum(Xsdpt3cvx.*U))-eps)^2);


%% Solving it with SDPT3
%
% Okays, so feeding it directly to SeDuMi makes it go a lot faster, thats
% good to confirm, but it seems that using SDPT3 is faster than SeDuMi

%[blk,At,C,bt,perm] = read_sedumi(Asedumi,bsedumi,csedumi,K);

blk{1,1}='s'; blk{1,2}= n;
blk{2,1}='s'; blk{2,2}= 2;
blk{3,1}='l'; blk{3,2}= 2;

bsdpt3=[b;eps;0;1];

for k=1:length(b)
    AA1{k}=mat(A(k,:));
    AA2{k}=[0,0;0,0];
end

AA1{length(b)+1}=-U;
AA1{length(b)+2}=zeros(n);
AA1{length(b)+3}=zeros(n);

AA2{length(b)+1}=[0, 0.5;0.5,0];
AA2{length(b)+2}=[0, 0.5;0.5,0];
AA2{length(b)+3}=[1,0;0,0];

At(1) = svec(blk(1,:),AA1);
At(2) = svec(blk(2,:),AA2);
At{3}= [zeros(length(b),2); -1,0;0,-1;0,0]';
C{1}=zeros(n);
C{2}=[0, mu/2;mu/2, c/2];
C{3}=[0,0];
[obj,X,y,Z] = sdpt3(blk,At,C,bsdpt3);

%%
% Big issue right now:  Why does Ax=0 instead of b?  This is really weird
% because im doing the same thing as the example code.

norm(At{1}'*svec(blk(1,:), X{1})+ At{2}'*svec(blk(2,:), X{2})+ At{3}'*svec(blk(3,:), X{3}))


%% Varying Noise
% Doing it correctly this time, last time i enforced equallity instead of
% saying ||Ax-b|| <= delta


clear;
if exist("Nov14Section1.mat", 'file')~=2
    rankTol=0.001;
    x=randn(10,50);
    D=x'*x;
    noise=0:0.03:0.15;
    for i=1:length(noise)
        
        [A,b] = sampleUniformSymmetric(D,0.55,noise(i));
        
        gamma=1.1*norm(A*vec(D)-b);
        
        [ Xnuc ] = nuclearNormPSD(zeros(50,50),A,b, gamma, [], [], []);
        [ Xar ] = aggregatedRelaxationPADM(zeros(50,50),A,b, gamma, [], [], [10, 0.95]);
        [ Xmpec ] = Rank_MPEC_PADM(zeros(50,50),A, b,sqrt(gamma),zeros(1,2500),[0],1e6);
        [ Xlogdet ] = logDet(zeros(50,50),A,b, gamma, [], [], [1,1]);
        [ Xmm, obj ]=matrixRelaxationMM(zeros(50,50),A,b, gamma, [], [], [1,1]);
        
        error(i,1)=norm(Xnuc-D, 'fro')/norm(D, 'fro');
        error(i,2)=norm(Xmpec-D, 'fro')/norm(D, 'fro');
        error(i,3)=norm(Xlogdet-D, 'fro')/norm(D, 'fro');
        error(i,4)=norm(Xar-D, 'fro')/norm(D, 'fro');
        error(i,5)=norm(Xmm-D, 'fro')/norm(D, 'fro');
        
        
        rank(i,1)=sum(svd(Xnuc)>rankTol);
        rank(i,2)=sum(svd(Xmpec)>rankTol);
        rank(i,3)=sum(svd(Xlogdet)>rankTol);
        rank(i,4)=sum(svd(Xar)>rankTol);
        rank(i,5)=sum(svd(Xmm)>rankTol);
    end
    save("Nov14Section1.mat")
else
    load("Nov14Section1.mat")
end

[rankTable, errorTable]=displayResults(noise, rank, error, [],{'Noise', 'Nuclear_Norm', 'MPEC_PADM', 'Log_Det', 'Aggregate_Relaxation', 'Matrix_Relaxation'});

%% PMU tests

if exist("nov14Section2.mat", 'file')~=2    
    plist=0.1:0.1:0.9;
    rank=[];
    error=[];
    
    load("pmuData.mat")
    data=precondition(data);
    dataReal=[real(data), imag(data)];
    rankTol=1;
    figure();
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [.1, 0.1, .9, 0.9])

    for i=1:length(plist)
        data=precondition(data);
        dataReal=[real(data), imag(data)];
        [ M,b,row,col ] = sampleUniform(data,plist(i));
        rowReal=[row, row];
        colReal=[col, col+37];
        bReal=[real(b), imag(b)];
        MReal=[real(M), imag(M)];
        [p,m]=size(data);
                

        
        tic;
        Xsvt=svt(rowReal, colReal, bReal, 600, 74, 100, .1, 0.0001, 1000);
        runtime(i,1)=toc;
        
        tic;
        [Xmr, k]=sIRLS(zeros(600, 74),MReal,-2,10000, 0.000001, 1000);
        runtime(i,2)=toc;
        error(i,1)=norm(Xsvt-dataReal, 'fro')/norm(dataReal, 'fro');
        error(i,2)=norm(Xmr-dataReal, 'fro')/norm(dataReal, 'fro');
        
        rank(i,1)=sum(svd(Xsvt)>rankTol);
        rank(i,2)=sum(svd(Xmr)>rankTol);
        
        subplot(ceil(length(plist)/3), 3, i)
        semilogy(svd(Xmr))
        hold on
        semilogy(svd(Xsvt))
        xlabel("i")
        ylabel("i^{th} Singular Value")
        title(plist(i)*100+"% data missing")
        legend(["SVT", "Matrix Relaxation"])
        
        
    end
    
    save("nov14Section2.mat")
else
    openfig('nov14section2.fig')
    load("nov14Section2.mat")
end
[rankTable, errorTable, runTimeTable]=displayResults(plist, rank, error, runtime, {'Percent_of_Data_Missing', 'SVT', 'Matrix_Relaxation'});

