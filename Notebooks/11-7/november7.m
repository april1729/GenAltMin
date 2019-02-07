%% Showing MM Algorithm is resiliant to choice of lambda.
clear;
if exist("Section1.mat", 'file')~=2
    d=3;
    n=30;
    rankTol=0.1;
    plist=0.06:0.03:0.3;
    for i=1:length(plist)
        [points, distances, A, b]=generateSensorProblem(n,d,plist(i), 0);
        
        [ Xmm1, obj ]=matrixRelaxationMM(zeros(n,n), [], [], 1, A, b, [0.01,1]);
        [ Xmm2, obj ]=matrixRelaxationMM(zeros(n,n), [], [], 1, A, b, [.1,1]);
        [ Xmm3, obj ]=matrixRelaxationMM(zeros(n,n), [], [], 1, A, b, [1,1]);
        [ Xmm4, obj ]=matrixRelaxationMM(zeros(n,n), [], [], 1, A, b, [10,0.95]);
        
        error(i,1)=calcErrorPoints(Xmm1, points);
        error(i,2)=calcErrorPoints(Xmm2, points);
        error(i,3)=calcErrorPoints(Xmm3, points);
        error(i,4)=calcErrorPoints(Xmm4, points);
        
        
        rank(i,1)=sum(svd(Xmm1)>rankTol);
        rank(i,2)=sum(svd(Xmm2)>rankTol);
        rank(i,3)=sum(svd(Xmm3)>rankTol);
        rank(i,4)=sum(svd(Xmm4)>rankTol);
    end
    save("section1.mat")
else
    load("section1.mat")
end
[rankTable, errorTable]=displayResults(plist, rank, error, {'Percent_of_Data_Known', 'a', 'b', 'c', 'd'});
%% Varying Noise
clear;
if exist("Section2.mat", 'file')~=2
    rankTol=0.001;
    x=randn(10,50);
    D=x'*x;
    noise=0:0.02:0.2;
    for i=1:length(noise)
        
        [A,b] = sampleUniformSymmetric(D,0.55,noise(i));
        
        
        [ Xnuc ] = nuclearNormPSD(zeros(50,50),[],[],0,A,b, []);
        [ Xar ] = aggregatedRelaxationPADM(zeros(50,50),[],[],0,A,b, [10, 0.95]);
        [ Xmpec ] = Rank_MPEC_PADM(zeros(50,50),zeros(1,2500),[0],1,A,b,1e6);
        [ Xlogdet ] = logDet(zeros(50,50),zeros(1,2500),[0],1,A,b,[]);
        [ Xmm, obj ]=matrixRelaxationMM(zeros(50,50), [], [], 1, A, b, [1,1]);
        
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
    save("section2.mat")
else
    load("section2.mat")
end

[rankTable, errorTable]=displayResults(noise, rank, error, {'Noise', 'Nuclear_Norm', 'MPEC_PADM', 'Log_Det', 'Aggregate_Relaxation', 'Matrix_Relaxation'});


%% Redoing experiement from last week with better parameters for Aggregatged Relaxation
if exist("Section3.mat", 'file')~=2
    
    rankTol=0.001;
    x=randn(10,50);
    D=x'*x;
    plist=[0.3:0.02:0.6];
    
    for i=1:length(plist)
        
        [A,b] = sampleUniformSymmetric(D,plist(i), 0);
        
        
        [ Xnuc ] = nuclearNormPSD(zeros(50,50),[],[],0,A,b, []);
        [ Xar ] = aggregatedRelaxationPADM(zeros(50,50),[],[],0,A,b, [10, 0.95]);
        [ Xmpec ] = Rank_MPEC_PADM(zeros(50,50),zeros(1,2500),[0],1,A,b,1e6);
        [ Xlogdet ] = logDet(zeros(50,50),zeros(1,2500),[0],1,A,b,[]);
        [ Xmm, obj ]=matrixRelaxationMM(zeros(50,50), [], [], 1, A, b, [1,1]);
        
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
    save("section3.mat")
else
    load("section3.mat")
end
[rankTable, errorTable]=displayResults(plist, rank, error,[], {'Percent_Data_Known', 'Nuclear_Norm', 'MPEC_PADM', 'Log_Det', 'Aggregate_Relaxation', 'Matrix_Relaxation'});

%% A slightly different formulation for Matrix Completion with PMU data
% Lets compare 2 methods: SVT and Matrix Relaxation, and lets try to
% seperate the real and imaginary parts of the entries, yeah, i like that.
%
% Here's exactly how im changing this matrix:

clear;
load("pmuData.mat")

newData=[real(data), imag(data)];

%%
% Heres the singular values of the matrix with and without seperating
% the real and imaginary parts
figure()
semilogy(svd(data))
hold on
semilogy(svd(newData));
xlabel("i")
ylabel("\sigma_i")
legend(["Old Matrix", "New Matrix"])
title("Singular value distribution for new and old matrix")

%% 
% To be more quantitative, let's show that both matricies are approximately
% rank 10.  This is important because now we have a 600x74 rank 10 matrix
% instead of a 600x37 rank 10 matrix.  The rank to column ratio is lower
% now.

[u,s,v]=svds([real(data), imag(data)], 10);
fprintf("Error in rank 10 approximation for seperating real and imag parts %f \n", norm([real(data), imag(data)]-u*s*v', 'fro')/norm([real(data), imag(data)], 'fro'))



[u,s,v]=svds([data], 10);
fprintf("Error in rank 10 approximation for original formulation %f \n", norm(data-u*s*v', 'fro')/norm(data, 'fro'))


%% Comparing Result
%
% Lets try to look at how well we can do MC on both of these matricies, I'm
% hoping it does better on this new fomulation.
%
% To make it a good example, what im going to do is sample the whole
% matrix, and then reconstruct the matrix based off of that so that im
% sampling the same stuff, and i'm not sampled the real but not imaginary
% part of a single entry.

if exist("Section4.mat", 'file')~=2
    data=precondition(data);
    dataReal=[real(data), imag(data)];
    [ M,b,row,col ] = sampleUniform( data,0.5 );
    rowReal=[row, row];
    colReal=[col, col+37];
    bReal=[real(b), imag(b)];
    [p,m]=size(data);
    X=svt(row, col, b,600, 37, 100, .1, 0.0001, 200);
    Xreal=svt(rowReal, colReal, bReal, 600, 74, 100, .1, 0.0001, 200);
    save("section4.mat")
else
    load("section4.mat")
end
fprintf("Error in new formulation: %f \n", norm(Xreal-dataReal, 'fro')/norm(dataReal, 'fro'))
fprintf("Error in old formulation: %f \n", norm([real(X), imag(X)]-dataReal, 'fro')/norm(dataReal, 'fro'))
%%
% Cool!! It works better!!! This makes a lot of sense.  I think this better
% takes advatage of the low dimensionality of the data.  I wonder how this
% effects the coherence
fprintf("Coherence of new formulation: %f \n", coherence(data))
fprintf("Coherence of old formulation: %f \n", coherence(dataReal))
%%
% These are both pretty awful, but I suppose this is a little bit better.