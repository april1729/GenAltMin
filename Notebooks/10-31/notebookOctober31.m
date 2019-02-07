%% October 24th-31st, 2018 Notebook
% My main goal this week is to get convincing computational results that this
% method works well. The algorithm is modified from last week in the sense
% that it now has epsilon decay over time, and so I started by testing that
% out.
%
% Since this seems to be performing well on the one case that ive given it,
% I want to try harder problems, so im going to include noise.  I'm also
% going to try doing the EDM promblem, and compare it to the LogDet and
%% Decrease Epsilon with each iteration
%


if exist("Section1_testDecay.mat", 'file')~=2
    pavg=0.70;
    eps=15;
    x=randn(10,50);
    D=x'*x;
    [ M,b,row,col ] = sampleUniform( D,pavg);
    P=~(M==0);
    [X1,U1, obj1]= aggregatedRelaxationPADM(M,P, eps, 1 );
    [X8,U8, obj8]= aggregatedRelaxationPADM(M,P, eps, 0.8);
    [X5,U5, obj5]= aggregatedRelaxationPADM(M,P, eps, 0.5);
    
    [u,s1,v]=svd(X1);
    [u,s8,v]=svd(X8);
    [u,s5,v]=svd(X5);
    save("Section1_testDecay.mat")
else
    load("Section1_testDecay.mat")
end
fprintf("Eps Decay Rate: %f \t Log Percent Error: %f \t log 11th Largest SV: %f \t Number of Iterations: %i\n", 1,log10(norm(D-X1, 'fro')/norm(D, 'fro')),log10(s1(11,11)), length(obj1));
fprintf("Eps Decay Rate: %f \t Log Percent Error: %f \t log 11th Largest SV: %f \t Number of Iterations: %i\n", 0.8 ,log10(norm(D-X8, 'fro')/norm(D, 'fro')),log10(s8(11,11)), length(obj8));
fprintf("Eps Decay Rate: %f \t Log Percent Error: %f \t log 11th Largest SV: %f \t Number of Iterations: %i\n", 0.5,log10(norm(D-X5, 'fro')/norm(D, 'fro')),log10(s5(11,11)), length(obj5));
%%
rankTol=0.001;

if exist("Section2.mat", 'file')~=2
    plist=0.3:0.01:0.6;
    for i=1:length(plist)
        pavg=plist(i);
        [M] = sampleUniformSymmetric(D,pavg);
        P=~(M==0);
        [ XNuc ] = nuclearNormPSD( M,P );
        relErrnuc(i)=norm(x'*x-XNuc,'fro')/norm(x'*x,'fro');
        rankResnuc(i)=sum(svd(XNuc)>rankTol);
        [X,U, obj]= aggregatedRelaxationPADM(M,P, eps,0.7);
        relErr(i)=norm(x'*x-X,'fro')/norm(x'*x,'fro');
        rankRes(i)=sum(svd(X)>rankTol);
        fprintf("Pavg: %f \t PADM Error: %f \t Nuc Error: %f \t PADM Iterations: %i \n", pavg,relErr(i), relErrnuc(i), length(obj));
    end
    save("Section2.mat")
else
    load("Section2.mat")
end
subplot(2,1,1)
plot(plist, rankRes)
hold on
plot(plist, rankResnuc)
legend(['Aggregate Relaxation', "Nuclear Norm with CVX"])
ylabel("Rank (0.001 tolerance)")

subplot(2,1,2)
plot(plist, relErr)
hold on
plot(plist, relErrnuc)
ylabel("Relative Error")
xlabel("Percent of data missing")
%%
% This looks to be working well, theres just one minor issue: If the number
% of iterations is to small, it doesnt get to a small enough value of
% epsilon to reduce the residual ion the singular calues past what it
% decides is okay to be big, wow im explaining this bad, but you know what
% i mean.  epsilon is too big so the 11th sv is too big.  yeah, that.  I
% just happened once here at the p=0.42.  Honestly i have no idea what
% happend at 0.4.  Not sure what a better convergence condition would be.
% Ask Prof Mitchell if he has any ideas for that!!  But the median of a lot
% of runs will look pretty good.

%% Adding noise!
% Whats a good way to add noise?  Well, um, it should be symmetric, right?
% Yeah, probably, because if this was a measurement, then id be measuring
% Mij and Mji to be the same thing, im not measuring two things, so theyll
% have the same noise.  Otherwise youre gaining information.
%
% One thing to think about, I'm going to need some kind of tolerance
% parameter, $ \gamma $ for the equality constraint
%
% $$||P_{\Omega}(X-M)||_F^2 \leq \gamma$$
%
% Whats a good choice for $\gamma$?  Well, it has to be larger than
% $||P_{\Omega}(N)||_F^2$ Where $N$ is the noise measurement.  I could do
% the math or just type that into matlab real quick:

noise=randn(50);
noise=0.05*.5*(noise'+noise);
disp(norm(P.*noise, 'fro')^2);

%%
% So then, maybe we could just say $\gamma=1$ for now?  Lets try that.
%
clear;
rankTol=0.5;
gamma=0.6;
eps=5;

x=randn(10,50);
noise=randn(50);
noise=.5*(noise'+noise);

D=x'*x+ 0.05*noise;

if exist("Section3.mat", 'file')~=2
    plist=0.4:0.1:0.7;
    for i=1:length(plist)
        pavg=plist(i);
        [M] = sampleUniformSymmetric(D,pavg);
        P=~(M==0);
        [ XNuc ] = nuclearNormPSD( M,P , gamma);
        relErrnuc(i)=norm(x'*x-XNuc,'fro')/norm(x'*x,'fro');
        rankResnuc(i)=sum(svd(XNuc)>rankTol);
        [X,U, obj]= aggregatedRelaxationPADM(M,P, 2.5,1, gamma);
        relErr(i)=norm(x'*x-X,'fro')/norm(x'*x,'fro');
        rankRes(i)=sum(svd(X)>rankTol);
        fprintf("Pavg: %f \t PADM Error: %f \t Nuc Error: %f \t PADM Iterations: %i \n", pavg,relErr(i), relErrnuc(i), length(obj));
    end
    save("Section3.mat")
else
    load("Section3.mat")
end

figure()
subplot(2,1,1)
plot(plist, rankRes)
hold on
plot(plist, rankResnuc)
legend(['Aggregate Relaxation', "Nuclear Norm with CVX"])
ylabel("Rank (0.5 tolerance)")

subplot(2,1,2)
plot(plist, relErr)
hold on
plot(plist, relErrnuc)
ylabel("Relative Error")
xlabel("Percent of data known")

%%
% After much trial and error, I made this work, however the method of
% decreasing epsilon at each iteration doesnt really work.  It just ends up
% converging to an higher objective function.  This indicates that because
% the case with noise takes longer to converge, we should let it converge
% and then increas the parameter.
%
% Idk, maybe ill do this sometime soon.  If so i'll put it right here!


%% Comparison to more methods
%
% So heres a good idea: lets adapt this code to instead solve this problem:
%
% $$ min_{X} rank(X)$$
%
% $$ ||Ax-b||_{2} \leq \delta $$
%
% $$ Ex = g $$
%
% $$ X \succeq 0 $$
%
% With this, i'll be able to go between euclidean distnace, matrix
% completion, or any other rank minimization problem.  All solvers should
% be called as solver(x0,A,b,delta,E,g, params), and params should be an
% optional input.
%
% Update: I did it.  So thats cool.  i also modified the code for the PADM
% in the non-relaxed version too.
clear;
if exist("Section4.mat", 'file')~=2
    
    rankTol=0.001
    x=randn(10,50);
    D=x'*x;
    plist=[0.3:0.01:0.6];
    
    for i=1:length(plist)
        
        [A,b] = sampleUniformSymmetric(D,plist(i));
        
        
        [ Xnuc ] = nuclearNormPSD(zeros(50,50),[],[],0,A,b, []);
        [ Xar ] = aggregatedRelaxationPADM(zeros(50,50),[],[],0,A,b, [10, 0.95]);
        [ Xmpec ] = Rank_MPEC_PADM(zeros(50,50),zeros(1,2500),[0],1,A,b,1e6);
        [ Xlogdet ] = LogDetHeuristic(zeros(50,50),zeros(1,2500),[0],1,A,b,1e6);
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
    save("section4.mat")
else
    load("section4.mat")
end


rankTable = array2table([plist', rank], 'VariableNames', {'Percent_of_Data_Known', 'Nuclear_Norm', 'MPEC_PADM', 'Log_Det', 'Augmented_Relaxation', 'Matrix_Relaxation'});
errorTable = array2table([plist', error], 'VariableNames', {'Percent_of_Data_Known', 'Nuclear_Norm', 'MPEC_PADM', 'Log_Det', 'Augmented_Relaxation', 'Matrix_Relaxation'});

disp("Rank:");
disp(rankTable);

disp("Error:");
disp(errorTable);


figure()
subplot(2,1,1)
plot(plist, rank)
legend(["Nuclear Norm", "MPEC PADM", "Log Det", "Augmented Relaxation", "Matrix Relaxation"],'Location','EastOutside')
ylabel("Rank (0.001 tolerance)")

subplot(2,1,2)
plot(plist, error)
ylabel("Relative Error")
xlabel("Percent of data missing")


%% Euclidean Distance Problem
%
% We now consider the Euclidean Distance Probelem.  Lets pose it like this:
% we know the location of one anchor point.  And that should be enough?
% no, no we need the location of d anchors where the points are in
% $R^d$.  Or, we could just say rotate it all after we do rank
% minimization?  Yeah, thats probably better?  I just need a way to compare
% these results is all.  so, then, we'll say this.  Solve the problem,
% rotate it so that it best fits where we think out anchor points are.
% nice, good plan.
clear;


if exist("Section5.mat", 'file')~=2
    
    d=3;
    n=30;
    rankTol=0.001;
    plist=0.06:0.03:0.3;
    
    for i=1:length(plist)
        [points, distances, A, b]=generateSensorProblem(n,d,plist(i), 0);
        
        [ Xnuc ] = nuclearNormPSD(zeros(n,n),[],[],0,A,b, []);
        [ Xar ] = aggregatedRelaxationPADM(zeros(n,n),[],[],0,A,b, [5, 0.9]);
        [ Xmpec ] = Rank_MPEC_PADM(zeros(n,n),zeros(1,n^2),[0],1,A,b,1e6);
        [ Xlogdet ] = LogDetHeuristic(zeros(n,n),zeros(1,n^2),[0],1,A,b,1e6);
        [ Xmm, obj ]=matrixRelaxationMM(zeros(n,n), [], [], 1, A, b, [1,1]);
        
        error(i,1)=calcErrorPoints(Xnuc, points);
        error(i,2)=calcErrorPoints(Xar, points);
        error(i,3)=calcErrorPoints(Xmpec, points);
        error(i,4)=calcErrorPoints(Xlogdet, points);
        error(i,5)=calcErrorPoints(Xmm, points);
        
        
        rank(i,1)=sum(svd(Xnuc)>rankTol);
        rank(i,2)=sum(svd(Xmpec)>rankTol);
        rank(i,3)=sum(svd(Xlogdet)>rankTol);
        rank(i,4)=sum(svd(Xar)>rankTol);
        rank(i,5)=sum(svd(Xmm)>rankTol);
    end
    
    save("section5.mat")
else
    load("section5.mat")
end
rankTable = array2table([plist', rank], 'VariableNames', {'Percent_of_Data_Known', 'Nuclear_Norm', 'MPEC_PADM', 'Log_Det', 'Augmented_Relaxation', 'Matrix_Relaxation'});
errorTable = array2table([plist', error], 'VariableNames', {'Percent_of_Data_Known', 'Nuclear_Norm', 'MPEC_PADM', 'Log_Det', 'Augmented_Relaxation', 'Matrix_Relaxation'});

disp("Rank:");
disp(rankTable);

disp("Error:");
disp(errorTable);


figure()
subplot(2,1,1)
plot(plist, rank)
legend(["Nuclear Norm", "MPEC PADM", "Log Det", "Augmented Relaxation", "Matrix Relaxation"], 'Location','EastOutside')
ylabel("Rank (0.001 tolerance)")

subplot(2,1,2)
plot(plist, error)
ylabel("Relative Error")
xlabel("Percent of data missing")


