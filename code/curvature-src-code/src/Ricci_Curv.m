function [all_curv] = Ricci_Curv(start_pt,end_pt,alpha)
%curv 
% alpha: the mass that remains at the same place for Markov random walk
% alpha>=0: fixed value alpha 
% alpha=-1: alpha=pi
% alpha=-2: alpha=1-pi
% alpha=-3: alpha=rho_i/sum_j~i(rho_j)

% -------------------------------------------------------

load('datasets_new/lusc_tcga.mat','adj','CNA','gene_list')
output_name = 'test';

p = gcp('nocreate');
if isempty(p)
    numcores = feature('numcores');
    myCluster = parcluster('local');
    myCluster.NumWorkers = numcores;
    saveProfile(myCluster); 
    parpool('local',numcores) 
end 

[Ng,Np]=size(CNA);% Ng: # of genes; Np: # of patients
[u,v]=find(triu(adj));
edge_list=sortrows([u,v],1);
Ne=size(edge_list,1);%Ne: # of edges

all_curv=zeros(Np,Ne);
all_curv_c=zeros(Np,Ne);
all_W = zeros(Np,Ne);

%parpool('local',40);
for pi=start_pt:end_pt
    rho=zeros(Ng,1);
    %inv_mea=zeros(Ng,Nc);
    rho=CNA(:,pi)+2.1; %remove 2.1 for RNA 
    if alpha==-1||alpha==-2
        inv_mea=invariant(adj,rho);
    else
        inv_mea=0;
    end
    Adj_w=zeros(size(adj));
    Adj_sw=zeros(size(adj));
    Adj_temp=adj.*rho(:)';
    sum_temp=adj*rho(:);
    Adj_temp=Adj_temp./sum_temp;
    Adj_w(:,:)=Adj_temp;
    Adj_temp=(Adj_temp+Adj_temp')/2;
    Adj_temp(Adj_temp~=0)=1./sqrt(Adj_temp(Adj_temp~=0));   
    Adj_sw(:,:)=Adj_temp;
    
    G=graph(Adj_sw);
    D=distances(G);
    parfor k=1:Ne
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
       	  alpha_i=rho(u)/sum(rho(setu));
       elseif alpha>=0
          alpha_i=alpha;
       end
       rho_u=[(1-alpha_i)*Adj_w(u,setu0),alpha_i];
       if alpha==-1
          alpha_i=inv_mea(v,1);
       elseif alpha==-2
          alpha_i=1-inv_mea(v,1);
       elseif alpha==-3
          alpha_i=rho(v)/sum(rho(setv));
       elseif alpha>=0
          alpha_i=alpha;
       end
       rho_v=[(1-alpha_i)*Adj_w(v,setv0),alpha_i];
       cc=Adj_sw(u,v);
       DD=D(setu,setv);
       d=Kantorovich(DD,rho_u,rho_v);
       curv=1-d/cc;
       all_curv(pi,k)=curv;
       all_curv_c(pi,k)=cc-d;
       all_W(pi,k)=d;
    end
    curv_i=all_curv(pi,:);
    curv_c_i=all_curv_c(pi,:);
    W_i=all_W(pi,:);
    pi
    save([output_name '_' num2str(alpha) '_' num2str(pi) '.mat'],'curv_i','curv_c_i','pi','W_i');
end
% save([output_name '_dynaAlpha_' num2str(gamma) '_all.mat'],'all_curv','all_curv_c','edge_list');
end

