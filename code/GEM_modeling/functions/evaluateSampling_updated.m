function [res] = evaluateSampling(model_ref, model, x, y, alpha, eff_thres, out_dir)
%EVALUATESAMPLING Summary of this function goes here
% Original code written by: Tilman Schäfers (2023)
% Updated by: Cyriel Huijer
% Last updated: 21 January 2025 
% Update: 
% - Fixed bug that computeCohen_d can only handle model prediction
%   smaller or equal than the number of reactions in the model. 
% - Log2FC included as extra output. 
% - Includes associated genes column
x=full(x);
y=full(y);
%Calculate set of overlapping reactions
rxn_overlap = intersect(model_ref.rxns, model.rxns);
disp("Total reactions ref/alt: " + numel(model_ref.rxns)+'/'+numel(model.rxns));
disp("Reaction overlap: " + numel(rxn_overlap));
x1 = x(getIndexes(model_ref, rxn_overlap, "rxns"),:);
y1 = y(getIndexes(model, rxn_overlap, "rxns"),:);
%Compare sample stats
[stats, pVals] = compareTwoSamplesStat(x1,y1);
disp(stats)
disp(pVals)
%Caluclate effect size
eff_sizes = cell(length(rxn_overlap),1);
for i=1:size(x1,1)
    eff_sizes{i,1} = abs(computeCohen_d(y1(i,:),x1(i,:)));
end
% select stats
idx = pVals.ks < alpha & abs(cell2mat(eff_sizes)) > eff_thres;
ks_pvals = pVals.ks(idx);
ks_stats = stats.ks(idx);
eff_filtered = eff_sizes(idx);
disp("Significantly altered reactions: " + numel(ks_pvals));
% Select significant reactions
rxn_changed = rxn_overlap(idx,:);
rxn_idx_changed = getIndexes(model_ref,rxn_changed,"rxns");
rxnNames_changed = model_ref.rxnNames(rxn_idx_changed);
subsystem_changed =  model_ref.subSystems(rxn_idx_changed);

% Identify mets_in and mets_out
mets_in = cell(numel(rxn_changed), 1);
mets_out = cell(numel(rxn_changed), 1);
associated_genes = cell(numel(rxn_changed), 1);

for i = 1:numel(rxn_idx_changed)
    rxn_idx = rxn_idx_changed(i);
    mets_in{i} = model_ref.mets(model_ref.S(:, rxn_idx) < 0); % Metabolites consumed
    mets_out{i} = model_ref.mets(model_ref.S(:, rxn_idx) > 0); % Metabolites produced
    associated_genes{i} = model_ref.grRules{rxn_idx};
end


disp(subsystem_changed);
% Create final result table
stats_x = calcSampleStats(x1(idx,:));
stats_y = calcSampleStats(y1(idx,:));
% Calculate log2fc:
log2fc = log2(stats_x.mean ./ stats_y.mean); % Log2 fold change
raw_diff = stats_x.mean - stats_y.mean;
% Write final results table
res = table(rxn_changed,rxnNames_changed, subsystem_changed,associated_genes, stats_x.mean, stats_y.mean, stats_x.std, stats_y.std, ks_stats, ks_pvals, eff_filtered,log2fc,raw_diff, mets_in,mets_out);
res.Properties.VariableNames = {'rxn','rxnName','subsystem','associated_genes','mean_sample_x','mean_sample_y','std_sample_x','std_sample_y','ks_stats','ks_pVal','cohen_effect','log2FC','raw_diff','mets_in','mets_out'};
res.Properties.VariableNames{'mean_sample_x'} = append('mean_sample_',model_ref.id);
res.Properties.VariableNames{'mean_sample_y'} = append('mean_sample_',model.id);
res.Properties.VariableNames{'std_sample_x'} = append('std_sample_',model_ref.id);
res.Properties.VariableNames{'std_sample_y'} = append('std_sample_',model.id);
res = sortrows(res,8,'ascend');
% return res and idx of changes reactions
writetable(res,fullfile(out_dir,strcat(model_ref.id,'_',model.id,'_','SamplingRxnsAltered.csv')),'Delimiter',',');  

end

