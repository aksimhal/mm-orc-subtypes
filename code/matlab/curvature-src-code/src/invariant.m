function  [W]  = invariant(Adj, Rho)
% Calculate Markov chain invariant measure
% Jiening Zhu, 07/11/2021

[m,n]=size(Rho);
W=zeros(size(Rho));
for i=1:n
    R=repmat(Rho(:,i),1,m);
    R_sum=sum(R.*Adj)'; 

    W(:,i)=Rho(:,i).*R_sum;
%     v=Rho(:,i).*R_sum;
%     v(v>mean(v)*10)=mean(v)*10;
%     W(:,i)=v;
end
W=W./sum(W);
end

