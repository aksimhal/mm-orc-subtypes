% Collate data
addpath '/multiple-myeloma/notebooks'

N_subjects = 669;

rna_data = load('[insert file name]', 'rna_data');
rna_data = rna_data.rna_data; 
rna_data = cell2mat(rna_data);

ADJ = load('[insert file name]', 'ADJ');
ADJ = ADJ.ADJ;
ADJ = ADJ.' | ADJ;


[u,v]=find(triu(ADJ));
edge_list=sortrows([u,v],1);

N_edges = length(edge_list); 
overall_curvature = zeros(N_edges, N_subjects);

output_name = '/curvature_';

for sample_n=1:N_subjects
    fn = [output_name num2str(sample_n) '_0_edges.mat'];
    load(fn);
    overall_curvature(:, sample_n)= all_curv(:);
    
end

save('[insert file name]', 'overall_curvature');