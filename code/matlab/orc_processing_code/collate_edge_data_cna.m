% Collate data

N_subjects = 669;

cna_data = load('[insert file name]', 'cna_data');
cna_data = cna_data.cna_data; 
cna_data = cell2mat(cna_data);

ADJ = load('[insert file name]', 'ADJ');
ADJ = ADJ.ADJ;
ADJ = ADJ.' | ADJ;


[u,v]=find(triu(ADJ));
edge_list=sortrows([u,v],1);

N_edges = length(edge_list); 
overall_curvature = zeros(N_edges, N_subjects);

output_name = '../data/output_cna/curvature_';

for sample_n=1:N_subjects
    fn = [output_name num2str(sample_n) '_0_edges.mat'];
    load(fn);
    overall_curvature(:, sample_n)= all_curv(:);
    
end

save('[insert file name]', 'overall_curvature');