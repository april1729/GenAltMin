m=100;
n=50;
r=4;
D=rand(m,n);
[U,S,V] = svds(D,r);
D=U*S*V';
numSampledList=0.05:0.05:0.8;
nucRank=[];
palmRank=[];
mmRank=[];
gdRank=[];
svtRank=[];
icmcRank=[];

pList=[];
palmError=[];
nucError=[];
svtError=[];
icmcError=[];
mmError=[];
gdError=[];


for pavg=numSampledList
    [ M,b,row,col ] = sampleUniform( D,pavg );
    [ X ] =  nuclearNormCVX( b, row, col, m,n);

    V=rand(m,n);
    bestRank=1000;
    bestH=2;
    
    X2=IRLS(M,M,-2,20,0.00001, 1000);
    Xgd=sIRLS(M,M,-2,20,0.00001, 1000);
    
    for g=flip(logspace(-1,1.2,20))
        X2=IRLS(X2,M,-2,g,0.00001, 1000);
        Xgd=sIRLS(Xgd,M,-2,g,0.00001, 1000);
    end

    
    
    for p=logspace(-3, 1, 400)
    [V,U , H, g, rank, pHist] = nonsymmetricPALM( V,row, col, b, p,100);
    if (sum(eig(V*V')>0.01)<bestRank) || ((sum(eig(V*V')>0.01)==bestRank) & (sum(sum(U.*(V'*V)))<bestH))
        bestRank=sum(svd(V)>0.1);
        bestp=p;
        bestError=norm(D-V,'fro')/norm(D,'fro');
        bestH=sum(sum(U.*(V'*V)));
        fprintf("p= %f \t final rank= %i \t n-trace(U)= %f \t <U, V'V>= %f \t Maximum percent error: %f \n",p,sum(eig(V*V')>0.01), n-trace(U), sum(sum(U.*(V'*V))), max(max(abs((V-D)./D))))

    end
    end
    nucRank=[nucRank,sum(svd(X)>0.1)];
    palmRank=[palmRank,bestRank];
    mmRank=[mmRank,sum(svd(X2)>0.1)];
    gdRank=[gdRank,sum(svd(Xgd)>0.1)];

    
    pList=[pList,bestp];
    palmError=[palmError, bestError];
    nucError=[nucError, norm(D-X,'fro')/norm(D,'fro')];
    mmError=[mmError, norm(D-X2,'fro')/norm(D,'fro')];
    gdError=[gdError, norm(D-Xgd,'fro')/norm(D,'fro')];
end


figure()
subplot(2,1,1)
plot(numSampledList, nucRank,numSampledList, palmRank,numSampledList, mmRank,numSampledList, gdRank)
xlabel("Percent of Points Known")
ylabel("Rank")
legend("Nuclear Norm solved with CVX",  "PALM", "MM", "Projected Gradient Descent")
subplot(2,1,2)
plot(numSampledList, nucError, numSampledList, palmError,numSampledList, min(1,mmError), numSampledList, gdError)
xlabel("Percent of Points Known")
ylabel("Relative Norm Error")