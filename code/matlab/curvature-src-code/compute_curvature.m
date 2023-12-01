function [all_curv, all_W] = compute_curvature(data_vector, adj, alpha)


[u,v]=find(triu(adj, 1));
edge_list=sortrows([u,v],1);
Ne=size(edge_list,1);%Ne: # of edges

if alpha==-1||alpha==-2
    inv_mea=invariant(adj,data_vector);
else
    inv_mea=0;
end

Adj_w=zeros(size(adj));
Adj_sw=zeros(size(adj));
Adj_temp=adj.*data_vector(:)';
sum_temp=adj*data_vector(:);
Adj_temp=Adj_temp./sum_temp;
Adj_w(:,:)=Adj_temp;
Adj_temp=(Adj_temp+Adj_temp')/2;
Adj_temp(Adj_temp~=0)=1./sqrt(Adj_temp(Adj_temp~=0));
Adj_sw(:,:)=Adj_temp;

G=graph(Adj_sw);
D=distances(G);
all_curv=zeros(Ne,1);
all_W=zeros(Ne,1);
for k=1:Ne
    u=edge_list(k,1);
    v=edge_list(k,2);
    setu0=find(adj(u,:)==1);
    setu=[setu0,u];
    setv0=find(adj(v,:)==1);
    setv=[setv0,v];
    if alpha==-1
        alpha_i=inv_mea(u,1);
    elseif alpha==-2
        alpha_i=1-inv_mea(u,1);
    elseif alpha==-3
        alpha_i=data_vector(u)/sum(data_vector(setu));
    elseif alpha>=0
        alpha_i=alpha;
    end
    rho_u=[(1-alpha_i)*Adj_w(u,setu0),alpha_i];
    if alpha==-1
        alpha_i=inv_mea(v,1);
    elseif alpha==-2
        alpha_i=1-inv_mea(v,1);
    elseif alpha==-3
        alpha_i=data_vector(v)/sum(data_vector(setv));
    elseif alpha>=0
        alpha_i=alpha;
    end
    rho_v=[(1-alpha_i)*Adj_w(v,setv0),alpha_i];
    cc=Adj_sw(u,v);
    DD=D(setu,setv);
    d=Kantorovich(DD,rho_u,rho_v);
    curv=1-d/cc;
    all_curv(k)=curv;
    all_W(k)=d;
end
end