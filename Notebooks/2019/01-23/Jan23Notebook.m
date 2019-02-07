%% January 23rd
% This week, i revisited the symmetric factorization I was working on a
% year ago.  It turns out, this is an extremely effective method.  I
% programmed it, iral.m%
% Also, as a note, ive updated the test problem generating code so theres
% no loops, this way it can generate very large probelems.  
n=100;
r=10;

[D, A, b]=generateMatrixCompletionProblemPSD(n,4,0.2, 0);
P0=rand(n,r);

lambda=1;
tic;[ X,obj ] = iral(P0, A,b,{@(x) lambda./(lambda+x), r,0.5,500,0.0001});
traceInvTime=toc;

tic;[ Xnuc,objnuc ] = iral(P0, A,b,{@(x) 1, 10,0.5,500,0.0001});
nucTime=toc;
fprintf("                           100X100, rank 4 with 20 percent entries\n")
fprintf("                        -----------------------------------------------\n")
fprintf("Nuclear Norm  || Error: %f \t Time: %f \t 5th Eigenvalue: %f\n", norm(Xnuc-D,'fro')/norm(D,'fro'), nucTime, min(svds(Xnuc,5)));
fprintf("Trace Inverse || Error: %f \t Time: %f \t 5th Eigenvalue: %f\n", norm(X-D,'fro')/norm(D,'fro'), traceInvTime, min(svds(X,5)));


%%
n=1000;
r=10;

[D, A, b]=generateMatrixCompletionProblemPSD(n,4,0.025, 0);
P0=rand(n,r);
    
lambda=1;
tic;[ X,obj ] = iral(P0, A,b,{@(x) lambda./(lambda+x), r,0.5,500,0.0001});
traceInvTime=toc;

tic;[ Xnuc,objnuc ] = iral(P0, A,b,{@(x) 1, 10,0.5,500,0.0001});
nucTime=toc;

svdX=svds(X,5);
svdXNuc=svds(Xnuc, 5);
fprintf("                         1000x1000, rank 4 with 2.5 percent entries\n")
fprintf("                        -----------------------------------------------\n")
fprintf("Nuclear Norm  || Error: %f \t Time: %f \t 5th Eigenvalue: %f\n", norm(Xnuc-D,'fro')/norm(D,'fro'), nucTime, min(svds(Xnuc,5)));
fprintf("Trace Inverse || Error: %f \t Time: %f \t 5th Eigenvalue: %f\n", norm(X-D,'fro')/norm(D,'fro'), traceInvTime, min(svds(X,5)));

%%
n=100;
d=3;
r=10;
lambda=1;
P0=rand(n,r);
[points, distances, A, b]=generateSensorProblem(n,d,0.5, 0);
tic;[ X,obj ] = iral(P0, A,b,{@(x) lambda./(lambda+x), r,0.5,500,0.0001});
traceInvTime=toc;
calcErrorPoints(X, points)


%%
load cow
Dist=calcDistance(pt*pt');
n=length(Dist);
d=3;
r=10;
lambda=1;
P0=rand(n,r);
[A, b] = sampleUniformGram(Dist,0.0075);


tic;[ X,obj ] = iral(P0, A,b,{@(x) lambda./(lambda+x), r,0.5,500,0.0001});
traceInvTime=toc;

tic;[ Xnuc,objnuc ] = iral(P0, A,b,{@(x) 1, r,0.5,500,0.0001});
nucTime=toc;

figure()
ptr=reconstructPoints(X, eye(3));
ViewMesh(ptr,trg);
title("Trace Inverse cow")

figure()
ptrXnuc=reconstructPoints(Xnuc, eye(3));
ViewMesh(ptrXnuc,trg);
title("Nuclear Norm Cow")
fprintf("                        Cow with 0.8 percent of distances \n")
fprintf("                        -----------------------------------------------\n")
fprintf("Nuclear Norm  || Error: %f \t Time: %f \t 4th Eigenvalue: %f\n", calcErrorPoints(Xnuc, pt), nucTime, min(svds(Xnuc,4)));
fprintf("Trace Inverse || Error: %f \t Time: %f \t 4th Eigenvalue: %f\n",calcErrorPoints(X, pt), traceInvTime, min(svds(X,4)));
%%
[F,pt]=stlread('sphere.stl');
Dist=calcDistance(pt*pt');
n=length(Dist);
d=3;
r=10;
lambda=1;
P0=rand(n,r);
[A, b] = sampleUniformGram(Dist,0.0025);


tic;[ X,obj ] = iral(P0, A,b,{@(x) lambda./(lambda+x), r,0.5,500,0.0001});
traceInvTime=toc;

tic;[ Xnuc,objnuc ] = iral(P0, A,b,{@(x) 1, r,0.5,500,0.0001});
nucTime=toc;

figure()
ptr=reconstructPoints(X, eye(3));
ViewMesh(ptr,F);
title("Trace Inverse cow")

figure()
ptrXnuc=reconstructPoints(Xnuc, eye(3));
ViewMesh(ptrXnuc,F);
title("Nuclear Norm Cow")
fprintf("                        Sphere with 0.25 percent of distances \n")
fprintf("                        -----------------------------------------------\n")
fprintf("Nuclear Norm  || Error: %f \t Time: %f \t 4th Eigenvalue: %f\n", calcErrorPoints(Xnuc, pt), nucTime, min(svds(Xnuc,4)));
fprintf("Trace Inverse || Error: %f \t Time: %f \t 4th Eigenvalue: %f\n",calcErrorPoints(X, pt), traceInvTime, min(svds(X,4)));
