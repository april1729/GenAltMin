function [X,Y,obj,time] = genAltMin_v2(M,Omega,opts)
[m,n]=size(M);

if isfield(opts, 'mu')
    mu=opts.mu;
else
    mu=1;
end
if isfield(opts, 'f')
    f=opts.f;
else
    f=@(x) 50^2./(50+x).^2;
end
if isfield(opts, 'maxIter')
    maxIter=opts.maxIter;
else
    maxIter=100;
end
if isfield(opts, 'xTol')
    xTol=opts.xTol;
else
    xTol=1e-5;
end
if isfield(opts, 'r')
    r=opts.r;
else
    r=10;
end

if isfield(opts, 'gamma')
    gamma=opts.gamma;
else
    gamma=100*norm(M,'fro')/(sqrt(sum(sum(Omega))/(m*n))*r);
end

if isfield(opts, 'fTol')
    fTol=opts.fTol;
else
    fTol=1e-6;
end

if isfield(opts, 'Zfull')
    Zfull=opts.Zfull;
else
    Zfull=0;
end

[row,col]=ind2sub([m,n], find(Omega==1));
Known=find(Omega==1);
data=M(Known);
X=randn(m,r);
Y=randn(n,r);
P=[X;Y];
W=eye(r);
alf = 0;  increment = 1;

Z = X*Y';  Res =  data - Z(Known);
Z(Known) = data;
S = sparse(row, col, data-partXY(X', Y', row, col, length(row))', m, n);
tic;
time=[];obj=[];
res=norm(Res);
rank_reset=0;
D=0;
for k=1:maxIter
    X0=X;Y0=Y;res0=res;Res0=Res;W0=W;
    if Zfull
        Z0=Z;
        X=(Z*Y)*pinv((W/mu+Y'*Y));
        X=(1+alf) * X - alf* X0;
        
        Y=(Z'*X)*(pinv(W/mu+X'*X));
        Y=(1+alf) * Y - alf* Y0;
        
        Z = X*Y';  Res =  data - Z(Known);
    else
        X=(X*(Y'*Y)+(1)*(S*Y))*pinv((W+Y'*Y));
        X=(1+alf) * X - alf * X0;
        
        Y=(Y*(X0'*X)+(1)*(S'*X))*(pinv(W+X'*X));
        Y=(1+alf)*Y - alf*Y0*(pinv(X'*X+W)*(X'*X0+W))';
        
        Res=data-partXY(X', Y', row, col, length(row))';
    end
    
    
    D0=D;
    [Vw,D]=eig(X'*X+Y'*Y);
    % (sum(1- gamma./(gamma+diag(D))))/(sum(1- gamma./(gamma+diag(D0))))
    %  if (sum(1- gamma./(gamma+diag(D))))/(sum(1- gamma./(gamma+diag(D0))))>0.98
    gamma0=gamma;
    gamma=max(0.75*gamma, 1e-4);
    % end
    time(k)=toc;
    obj(k)=1000*(sum(1- 10./(10+diag(D))))+0.5*norm(Res)^2;
    res = norm(Res);
    if k>1; ratio = obj(k)/obj(k-1); else; ratio=0; end
    
    % adjust alf
    if rank_reset; rank_reset=0;
    else
        if ratio > 1
            increment = max([0.1*alf, 0.1*increment,0.5]);
            X = X0; Y = Y0; Res = Res0; res = res0; obj(k)=obj(k-1);
            alf = 0;    if Zfull; Z = Z0;end
        elseif ratio > 0.7
            increment = max(increment, 0.25*alf);
            alf = alf+ increment;
        end
    end
%     if r>sum(diag(D)>0.01)
%         r=sum(diag(D)>0.01)
%         %X0=X;Y0=Y;
%         [Ux,Sx,Vx]=svds(X,r);
%         [Uy,Sy,Vy]=svds(Y,r);
%         middle_matrix=(Sx*Vx'*Vy*Sy)^0.5;
%         X=Ux*middle_matrix;
%         Y=Uy*middle_matrix';
%         %norm(X*Y'-X0*Y0')
%         rank_reset=1;
%         alf=0;
%         increment=1;
%         Res=data-partXY(X', Y', row, col, length(row))';
%     end
    %rank_reset
    
    if Zfull; Z(Known) = data; else updateSval(S, Res, length(data)); end
    if(k <= 4)
        delta = inf;
    else
        delta = obj(k - 4:k - 1) - obj(k - 3:k);
        delta = abs(delta)/obj(k);
        delta = mean(delta)/1.5;
        % delta = max(delta);
    end
    
    
    [Vw,D]=eig(X'*X+Y'*Y);
    
    W=10*gamma*Vw*diag((gamma./(gamma+diag(D)).^2))*Vw';
    if rank_reset; W=10*gamma*eye(r);end
    %W=(alf +1)*W-alf*W0;

    %W=0;
    
    fprintf('it: %5i rk: %4d, rel. %3.1e r. %4.4f chg: %3.1e alf: %3.1e gamma: %3.1e delta: %3.1e\n ',...
        k,sum(diag(D)>0.01), res,ratio,obj(k),alf,gamma, delta);
    
    if delta<1e-8
        break
    end
    
    
    
    
end

end
