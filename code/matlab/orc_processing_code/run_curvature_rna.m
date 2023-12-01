function run_curvature_rna(start_subject, end_subject)

rna_data = load('[insert file name]', 'rna_data');
rna_data = rna_data.rna_data; 
rna_data = cell2mat(rna_data);

ADJ = load('[insert file name]', 'ADJ');
ADJ = ADJ.ADJ;
ADJ = ADJ.' | ADJ;


addpath(genpath('../'));
alpha = 0;

parfor n = start_subject:end_subject
    data_vector = rna_data(:, n) + 0.001;
    
    [all_curv, all_W] = compute_curvature(data_vector, ADJ, alpha);
    parsave(['output_rna/curvature' '_' num2str(n) '_' num2str(0) '_edges.mat'], all_curv, all_W);
end
end


function parsave(fname, all_curv, all_W)
save(fname, 'all_curv','all_W');
end
