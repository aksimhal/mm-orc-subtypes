function [d,PI] = Kantorovich(C,p0,p1)

[m,n]=size(C);
% %% cvx
% cvx_begin quiet
%     variable PI(m,n) nonnegative
%     minimize dot(C(:),PI(:))
%     subject to
%         PI*ones(n,1) == p0';
%         transpose(PI)*ones(m,1) == p1';
% cvx_end
% d=cvx_optval;

%% linprog
list1=1:m:(1+m*(n-1));
list2=1:m;
zerorow=zeros(1,m*n);
A=[];
for i=1:m
    row_temp=zerorow;
    row_temp(list1+(i-1))=1;
    A=[A;row_temp];
end
for i=1:n
    row_temp=zerorow;
    row_temp(list2+(i-1)*m)=1;
    A=[A;row_temp];    
end
b=[p0';p1'];
options = optimoptions('linprog','Display','none');
[PI,d]=linprog(C(:),[],[],A,b,zeros(m*n,1),[],options);
end

