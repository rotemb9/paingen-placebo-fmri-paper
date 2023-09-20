% This script reads the data from paingen placebo study (demporaphic, behavioral and
% brain measures), and then runs mixed effects models for the placebo and stimulus intensity effects, accounting for the
% familial structure of the data (consisting of monozygotic and dizygotic twins).

%% read the data - CHANGE DIRS TO YOUR LOCAL PATHS
main_dir = pwd; % change the main_dir to the local path where you put the data
z_score = true; % if true, uses behavioral and brain measures that were z scored within modality, to normalized (beta) estimates in the model. Will not affect the p value.

if z_score
    data = readtable(fullfile(main_dir, 'pain_evoked_brain_measures_z.csv'));
else
    data = readtable(fullfile(main_dir, 'pain_evoked_brain_measures.csv'));    
end

%% define outcomes to run models on
% running on all behavioral and brain measures
ind_gray_matter = find(strcmp(data.Properties.VariableNames, 'gray_matter'));
brain_vars = data.Properties.VariableNames(ind_gray_matter:end);
outcome = [{'Yint', 'Yunp'}, brain_vars]; % add pai intensity and pain unpleasantness measures to the brain measures

%% run the models
[mdl_heat, STATS_heat, mdl_press, STATS_press, mdl_heat_int, STATS_heat_int, mdl_press_int, STATS_press_int] = deal(cell(1,length(outcome)));
for i = 1:length(outcome)   
    disp(outcome{i});
    % thermal pain rating
    mdl_heat{i} = fitlme(data(data.heat > 0,:), sprintf(['%s ~ stimLvl + prodicaine + (stimLvl + prodicaine | family_ID) + ',...
        '(stim_mz + prodicaine_mz + int_mz -1 | participant_ID) + ',...
        '(stim_dz + prodicaine_dz + int_dz -1 | participant_ID)'], outcome{i}),'FitMethod','REML');
    [~,~,STATS_heat{i}] = fixedEffects(mdl_heat{i},'dfmethod','Satterthwaite');
    % mechanical pain rating
    mdl_press{i} = fitlme(data(data.heat <= 0,:), sprintf(['%s ~ stimLvl + prodicaine + (stimLvl + prodicaine | family_ID) + ',...
        '(stim_mz + prodicaine_mz + int_mz -1 | participant_ID) + ',...
        '(stim_dz + prodicaine_dz + int_dz -1 | participant_ID)'], outcome{i}),'FitMethod','REML');
    [~,~,STATS_press{i}] = fixedEffects(mdl_press{i},'dfmethod','Satterthwaite');    
    
    % same models with interaction
     mdl_heat_int{i} = fitlme(data(data.heat > 0,:), sprintf(['%s ~ stimLvl * prodicaine + (stimLvl + prodicaine | family_ID) + ',...
        '(stim_mz + prodicaine_mz + int_mz -1 | participant_ID) + ',...
        '(stim_dz + prodicaine_dz + int_dz -1 | participant_ID)'], outcome{i}),'FitMethod','REML');

    [~,~,STATS_heat_int{i}] = fixedEffects(mdl_heat_int{i},'dfmethod','satterthwaite');

    mdl_press_int{i} = fitlme(data(data.heat <= 0,:), sprintf(['%s ~ stimLvl * prodicaine + (stimLvl + prodicaine | family_ID) + ',...
        '(stim_mz + prodicaine_mz + int_mz -1 | participant_ID) + ',...
        '(stim_dz + prodicaine_dz + int_dz -1 | participant_ID)'], outcome{i}),'FitMethod','REML');

    [~,~,STATS_press_int{i}] = fixedEffects(mdl_press_int{i},'dfmethod','satterthwaite');           

end

%% create a table with stats
stats_table_condition = table;
stats_table_condition.outcome = repmat(outcome',2,1);
stats_table_condition.modality = repelem({'thermal';'mechanical'},length(outcome));
for ind = 1:length(outcome)
    % placebo
    stats_table_condition.estimate_prodicaine(ind) = round(STATS_heat{ind}.Estimate(3),3);
    stats_table_condition.estimate_prodicaine(ind + length(outcome)) = round(STATS_press{ind}.Estimate(3),3);
    stats_table_condition.SE_prodicaine(ind) = round(STATS_heat{ind}.SE(3),3);
    stats_table_condition.SE_prodicaine(ind + length(outcome)) = round(STATS_press{ind}.SE(3),3);
    stats_table_condition.t_prodicaine(ind) = round(STATS_heat{ind}.tStat(3),3);
    stats_table_condition.t_prodicaine(ind + length(outcome)) = round(STATS_press{ind}.tStat(3),3);
    stats_table_condition.DF_prodicaine(ind) = round(STATS_heat{ind}.DF(3),2);
    stats_table_condition.DF_prodicaine(ind + length(outcome)) = round(STATS_press{ind}.DF(3),2);
    stats_table_condition.p_prodicaine(ind) = round(STATS_heat{ind}.pValue(3),4);
    stats_table_condition.p_prodicaine(ind + length(outcome)) = round(STATS_press{ind}.pValue(3),4);
    stats_table_condition.lower_prodicaine(ind) = round(STATS_heat{ind}.Lower(3),4);
    stats_table_condition.lower_prodicaine(ind + length(outcome)) = round(STATS_press{ind}.Lower(3),4);
    stats_table_condition.upper_prodicaine(ind) = round(STATS_heat{ind}.Upper(3),4);
    stats_table_condition.upper_prodicaine(ind + length(outcome)) = round(STATS_press{ind}.Upper(3),4);

    % stimLvl
    stats_table_condition.estimate_stimLvl(ind) = round(STATS_heat{ind}.Estimate(2),3);
    stats_table_condition.estimate_stimLvl(ind + length(outcome)) = round(STATS_press{ind}.Estimate(2),3);
    stats_table_condition.SE_stimLvl(ind) = round(STATS_heat{ind}.SE(2),3);
    stats_table_condition.SE_stimLvl(ind + length(outcome)) = round(STATS_press{ind}.SE(2),3);
    stats_table_condition.t_stimLvl(ind) = round(STATS_heat{ind}.tStat(2),3);
    stats_table_condition.t_stimLvl(ind + length(outcome)) = round(STATS_press{ind}.tStat(2),3);
    stats_table_condition.DF_stimLvl(ind) = round(STATS_heat{ind}.DF(2),2);
    stats_table_condition.DF_stimLvl(ind + length(outcome)) = round(STATS_press{ind}.DF(2),2);
    stats_table_condition.p_stimLvl(ind) = round(STATS_heat{ind}.pValue(2),4);
    stats_table_condition.p_stimLvl(ind + length(outcome)) = round(STATS_press{ind}.pValue(2),4);
    stats_table_condition.lower_stimLvl(ind) = round(STATS_heat{ind}.Lower(2),4);
    stats_table_condition.lower_stimLvl(ind + length(outcome)) = round(STATS_press{ind}.Lower(2),4); 
    stats_table_condition.upper_stimLvl(ind) = round(STATS_heat{ind}.Upper(2),4);
    stats_table_condition.upper_stimLvl(ind + length(outcome)) = round(STATS_press{ind}.Upper(2),4); 
    
    % interaction
    stats_table_condition.estimate_int(ind) = round(STATS_heat_int{ind}.Estimate(end),3);
    stats_table_condition.estimate_int(ind + length(outcome)) = round(STATS_press_int{ind}.Estimate(end),3);
    stats_table_condition.SE_int(ind) = round(STATS_heat_int{ind}.SE(end),3);
    stats_table_condition.SE_int(ind + length(outcome)) = round(STATS_press_int{ind}.SE(end),3);
    stats_table_condition.t_int(ind) = round(STATS_heat_int{ind}.tStat(end),3);
    stats_table_condition.t_int(ind + length(outcome)) = round(STATS_press_int{ind}.tStat(end),3);
    stats_table_condition.DF_int(ind) = round(STATS_heat_int{ind}.DF(end),2);
    stats_table_condition.DF_int(ind + length(outcome)) = round(STATS_press_int{ind}.DF(end),2);
    stats_table_condition.p_int(ind) = round(STATS_heat_int{ind}.pValue(end),4);
    stats_table_condition.p_int(ind + length(outcome)) = round(STATS_press_int{ind}.pValue(end),4);
    stats_table_condition.lower_int(ind) = round(STATS_heat_int{ind}.Lower(end),4);
    stats_table_condition.lower_int(ind + length(outcome)) = round(STATS_press_int{ind}.Lower(end),4); 
    stats_table_condition.upper_int(ind) = round(STATS_heat_int{ind}.Upper(end),4);
    stats_table_condition.upper_int(ind + length(outcome)) = round(STATS_press_int{ind}.Upper(end),4); 
end

% placebo
stats_table_condition.direction_prodicaine(stats_table_condition.estimate_prodicaine > 0) = {'+'};
stats_table_condition.direction_prodicaine(stats_table_condition.estimate_prodicaine < 0) = {'-'};
stats_table_condition.significance_prodicaine(stats_table_condition.p_prodicaine >= 0.05) = {'n.s.'};
stats_table_condition.significance_prodicaine(stats_table_condition.p_prodicaine < 0.05) = {'*'};
stats_table_condition.significance_prodicaine(stats_table_condition.p_prodicaine < 0.01) = {'**'};
stats_table_condition.significance_prodicaine(stats_table_condition.p_prodicaine < 0.001) = {'***'};

% stim level
stats_table_condition.direction_stimLvl(stats_table_condition.estimate_stimLvl > 0) = {'+'};
stats_table_condition.direction_stimLvl(stats_table_condition.estimate_stimLvl < 0) = {'-'};
stats_table_condition.significance_stimLvl(stats_table_condition.p_stimLvl >= 0.05) = {'n.s.'};
stats_table_condition.significance_stimLvl(stats_table_condition.p_stimLvl < 0.05) = {'*'};
stats_table_condition.significance_stimLvl(stats_table_condition.p_stimLvl < 0.01) = {'**'};
stats_table_condition.significance_stimLvl(stats_table_condition.p_stimLvl < 0.001) = {'***'};

% int
stats_table_condition.significance_int(stats_table_condition.p_int >= 0.05) = {'n.s.'};
stats_table_condition.significance_int(stats_table_condition.p_int < 0.05) = {'*'};
stats_table_condition.significance_int(stats_table_condition.p_int < 0.01) = {'**'};
stats_table_condition.significance_int(stats_table_condition.p_int < 0.001) = {'***'};

if z_score
    writetable(stats_table_condition, fullfile(main_dir, 'stats_condition_models_z.csv'));
else
    writetable(stats_table_condition, fullfile(main_dir, 'stats_condition_models.csv'));
end

%% save models to a mat file
if z_score
    save(fullfile(main_dir, 'models_z'), 'mdl_heat', 'mdl_press', 'STATS_heat','STATS_press', 'mdl_heat_int', 'STATS_heat_int', 'mdl_press_int', 'STATS_press_int', 'outcome')
else
    save(fullfile(main_dir, 'models'), 'mdl_heat', 'mdl_press', 'STATS_heat','STATS_press', 'mdl_heat_int', 'STATS_heat_int', 'mdl_press_int', 'STATS_press_int', 'outcome')
end

%% organize table for paper
% placebo
table_paper_condition_placebo = table;
table_paper_condition_placebo.outcome = repelem(outcome',2);
table_paper_condition_placebo.modality = repmat({'thermal';'mechanical'},length(outcome),1);
stats_table_condition_matched = join(table_paper_condition_placebo, stats_table_condition,'Keys',{'outcome', 'modality'});
table_paper_condition_placebo.estimate = round(stats_table_condition_matched.estimate_prodicaine,3);
table_paper_condition_placebo.CI = strcat('[', num2str(round(stats_table_condition_matched.lower_prodicaine,3)), ', ', num2str(round(stats_table_condition_matched.upper_prodicaine,3)), ']');
table_paper_condition_placebo.SE = round(stats_table_condition_matched.SE_prodicaine,3);
table_paper_condition_placebo.t = round(stats_table_condition_matched.t_prodicaine, 2);
table_paper_condition_placebo.df = round(stats_table_condition_matched.DF_prodicaine,1);
table_paper_condition_placebo.p = round(stats_table_condition_matched.p_prodicaine,3);

% stim level
table_paper_condition_stimLvl = table;
table_paper_condition_stimLvl.outcome = repelem(outcome',2);
table_paper_condition_stimLvl.modality = repmat({'thermal';'mechanical'},length(outcome),1);
table_paper_condition_stimLvl.estimate = round(stats_table_condition_matched.estimate_stimLvl,3);
table_paper_condition_stimLvl.CI = strcat('[', num2str(round(stats_table_condition_matched.lower_stimLvl,3)), ', ', num2str(round(stats_table_condition_matched.upper_stimLvl,3)), ']');
table_paper_condition_stimLvl.SE = round(stats_table_condition_matched.SE_stimLvl,3);
table_paper_condition_stimLvl.t = round(stats_table_condition_matched.t_stimLvl, 2);
table_paper_condition_stimLvl.df = round(stats_table_condition_matched.DF_stimLvl,1);
table_paper_condition_stimLvl.p = round(stats_table_condition_matched.p_stimLvl,3);

% combine the effects to a single table
table_paper_condition_stimLvl.effect(:) = {'Stim level'};
table_paper_condition_placebo.effect(:) = {'Placebo'};
table_paper_condition_stimLvl.CI = cellstr(table_paper_condition_stimLvl.CI);
table_paper_condition_placebo.CI = cellstr(table_paper_condition_placebo.CI);
table_paper_condition = [table_paper_condition_stimLvl; table_paper_condition_placebo];
table_paper_condition = sortrows(table_paper_condition, {'outcome', 'modality', 'effect'}, {'ascend', 'descend', 'descend'});
table_paper_condition.modality(strcmp(table_paper_condition.modality, 'pressure')) = {'mechanical'};
if z_score
    writetable(table_paper_condition, fullfile(main_dir, 'table_paper_condition_z.csv'));
else
    writetable(table_paper_condition, fullfile(main_dir, 'table_paper_condition.csv'));
end