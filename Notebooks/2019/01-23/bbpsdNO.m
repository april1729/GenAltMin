function [ P , obj] = bbpsd(P0,W,A,b,params )
%BBALPSD Summary of this function goes here
%   Detailed explanation goes here
%   This solves the problem:
%
%   min trace(W*P'*P)+mu/2 ||A*vec(P*P')-b||^2
[n,r]=size(P0);
mu=params{3};
maxIter=params{4};
eps=params{5};
[~,j] = find(A);
[rows,cols]=ind2sub([n,n], j);
%P=P0;

%f=@(P) abs(r-trace(pinv(P'*P/lambda+eye(r)))+(mu/2)*norm(A*vec(P*P')-b,'fro')^2);
%f=@(P) trace(W*P'*P)+(mu/2)*norm(A*vec(P*P')-b,'fro')^2;%+10*norm(P-P0,'fro')^2;
%g=@(P) 2*P+2*mu*reshape(A'*(A*vec(P*P')-b),[n,n] )*P;%+20*(P-P0);


options = optimoptions('fminunc','SpecifyObjectiveGradient',true);
P= fminunc(@fun,P0,options);

% for k=1:maxIter
%     sold=s;
%     s=2*P*W+2*mu*reshape(A'*(A*vec(P*P')-b),[n,n])*P;
%     %s=2/lambda*P*inv(P'*P/lambda+eye(r))^2+2*mu*reshape(A'*(A*vec(P*P')-b),[n,n] )*P;
%     if k==1
%         alpha=10e-5;
%     else
%         %alpha=0.001*(vec(P-Pold)'*vec(P-Pold))/(vec(P-Pold)'*vec(s-sold));
%         alpha=10e-5;
%     end
%     Pold=P;
%
%     P=P-alpha*s;
%
%     obj(k)=f(P);
%
%     if k>2
%         norm((P-Pold), 'fro')/norm(Pold, 'fro');
%         if norm((P-Pold), 'fro')/norm(Pold, 'fro')<eps
%             break
%         end
%     end
%
% end
X=P*P';



    function [f,g]=fun(P)
        
        f= trace(W*P'*P)+(mu/2)*norm(A*vec(calculateSpaseX(P,rows,cols))-b,'fro')^2;%+10*norm(P-P0,'fro')^2;
        g= 2*P+2*mu*reshape(A'*(A*vec(calculateSpaseX(P,rows,cols))-b),[n,n] )*P;%+20*(P-P0);
    end

end
