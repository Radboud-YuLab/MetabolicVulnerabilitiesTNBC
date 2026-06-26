%% script for random sampling
% written by: Cyriel Huijer
% last updated: 16 Dec 2025


initCobraToolbox;
changeCobraSolver('gurobi');
setRavenSolver('gurobi');

clear
clc

cd C:\Users\chuijer\Documents\nextcloud_backup\PhD\Projects\Project1\1_manuscript\code\

load("data/modeling/MCF10A_model.mat")
MCF10A_model = EV_0_thr1_tINIT_model;
% measurements are in mmol/gDW
lb = readtable("data/exchange_fluxes/exch_fluxes_3D/for_modeling/20251217_lb.csv");
ub = readtable("data/exchange_fluxes/exch_fluxes_3D/for_modeling/20251217_ub.csv");

%% comment in to get % error runs
% % Extract numeric parts (all columns except 'metabolite')
% lb_vals = lb{:, 2:end};
% ub_vals = ub{:, 2:end};
% 
% % Calculate mean
% mean_vals = (lb_vals + ub_vals) / 2;
% 
% % Create output table
% mean_fluxes = array2table(mean_vals, ...
%     'VariableNames', lb.Properties.VariableNames(2:end));
% 
% % Add metabolite column back
% mean_fluxes = addvars(mean_fluxes, lb.metabolite, ...
%     'Before', 1, 'NewVariableNames', 'metabolite');
% numeric_data = mean_fluxes{:, 2:end};
% lb = mean_fluxes;
% ub = mean_fluxes;
% lb{:, 2:end} = numeric_data - abs((numeric_data * 0.2));
% ub{:, 2:end} = numeric_data + abs((numeric_data * 0.2));

conditions = lb.Properties.VariableNames(2:end) 
n_conditions = width(lb) - 1;
models = repmat({MCF10A_model}, n_conditions, 1);  % 18x1 cell array of models

for i = 1:length(models)
    % Display the current condition
    disp(conditions{i});
    % Apply the setHamsMedium_CH function each model
    models{i,1} = setHamsMedium_CH(models{i,1});
    % Find the biomass reaction index
    biomass_rxn_idx = find(strcmp(models{i,1}.rxns, 'MAR13082'));
    % Set the biomass objective coefficient to 1
    models{i,1}.c(biomass_rxn_idx) = 1;
    % Solve the LP problem for this model
    sol = solveLP(models{i,1}, 'max');
    % Display the optimal objective value
    disp(sol.f);
end

% Find rows where metabolite is methionine and valine, these were also
% removed in REGP manuscript as they are infeasible
rowsToRemove = ismember(lb.metabolite, {'methionine','valine'});

% Remove those rows
lb(rowsToRemove, :) = [];
ub(rowsToRemove, :) = [];

% optimize biomass rxn
max_solutions = [];
min_solutions = [];
for i = 1:length(conditions)
    cond = conditions{i};
    disp(cond);
    % Set exchange bounds for each condition
    models{i,1} = setExchangeBounds(models{i,1}, ...
                                  lb.metabolite, ...
                                  lb.(cond), ...
                                  ub.(cond), ...
                                  false);
    sol_max = optimizeCbModel(models{i},'max');
    sol_min = optimizeCbModel(models{i},'min');
    disp(sol_max.f);
    disp(sol_min.f)
    % if the solution is infeasible
    if sol_max.stat == 0
        max_solutions = [max_solutions,sol_max.f];
    end
    if sol_max.f >= 0
        max_solutions = [max_solutions,sol_max.f];
    end
    if sol_min.stat == 0
        min_solutions = [min_solutions,sol_min.f];
        disp("Infeasible");
    end
    if sol_min.f >= 0
        min_solutions = [min_solutions,sol_min.f];
    end
    models{i}.id = cond;
end

% plot max gr

figure;
bar(max_solutions);

% Set x-axis labels
xticks(1:numel(conditions));
xticklabels(conditions);

% Improve readability
xtickangle(45);

ylabel('Growth rate');
title('Max growth rates per condition');

%save("output/20251217_models_constrained.mat","models");

clear sol_max sol_min sol cond EV_0_thr1_tINIT_model i max_solutions min_solutions rowsToRemove 

%% constrain with GR
gr_table = readtable("data/growth_rate/celltiterglo_3D/output/20251217_growth_rate.csv");

for i = 1:length(conditions)
    models{i,1} = setParam(models{i,1},'lb',biomass_rxn_idx,gr_table{i,"lb"});
    models{i,1} = setParam(models{i,1},'ub',biomass_rxn_idx,gr_table{i,"ub"});
    sol = optimizeCbModel(models{i,1},'max');
    disp(sol.f);
end

%save("20251217_models_constrained.mat","models");

%% Run random sampling

% Run rs 1x so good_rxns is formed, this is possible since the same tINIT
% GEM is used. 
[~,good_rxns] = randomSampling(models{10},1,true, false, true, [], false);
solutions = cell(n_conditions, 1);   % preallocate cell array

for i = 1:nModels
    fprintf('Sampling model %d / %d: %s\n',i, nModels, conditions{i});
    fname = sprintf('data/modeling/rs_results/rs_solution_%02d_%s.mat', ...
                    i, conditions{i});
    solutions{i} = randomSampling(models{i},1000, true, false,true, good_rxns, false);
    sol_i = solutions{i};
    save(fname,'sol_i')
end

