% This script reads the data from paingen placebo study (demporaphic, behavioral and
% brain measures), and then runs mixed effects models for the correlations between behavioral and neural placebo-induced reductions, accounting for the
% familial structure of the data (consisting of monozygotic and dizygotic twins).
% written by Rotem Botvinik-Nezer

%% read the data
main_dir = pwd; % change the main_dir to the local path where you put the data
z_score = true; % if true, uses behavioral and brain measures that were z scored within modality, to normalized (beta) estimates in the model. Will not affect the p value.

if z_score
    analgesia_table = readtable(fullfile(main_dir, 'pain_evoked_analgesia_z.csv'));
else
    analgesia_table = readtable(fullfile(main_dir, 'pain_evoked_analgesia.csv'));    
end


%% define outcomes to run models on
% running on all brain measures
ind_gray_matter = find(strcmp(analgesia_table.Properties.VariableNames, 'gray_matter'));
outcome_analgesia = analgesia_table.Properties.VariableNames(ind_gray_matter:end);
outcome_analgesia = outcome_analgesia(~strcmp(outcome_analgesia, 'rank_Yint'));

% create variables for random intercepts and slopes
analgesia_table.zygosity = (contains(analgesia_table.twin_ID, 'DZ') | contains(analgesia_table.twin_ID, 'OS')) + 1; % zygosity coding: 1=monozygotic, 2=dizygotic
analgesia_table.Yint_mz(analgesia_table.zygosity == 1) = analgesia_table.Yint(analgesia_table.zygosity == 1);
analgesia_table.Yint_dz(analgesia_table.zygosity == 2) = analgesia_table.Yint(analgesia_table.zygosity == 2);
analgesia_table.rank_Yint_mz(analgesia_table.zygosity == 1) = analgesia_table.rank_Yint(analgesia_table.zygosity == 1);
analgesia_table.rank_Yint_dz(analgesia_table.zygosity == 2) = analgesia_table.rank_Yint(analgesia_table.zygosity == 2);

%% run the models
[mdl_heat_analgesia, STATS_heat_analgesia, mdl_press_analgesia, STATS_press_analgesia] = deal(cell(1,length(outcome_analgesia)));
if z_score
    warning('note that the rank measures are not z scored');
end
for i = 1:length(outcome_analgesia)
    disp(outcome_analgesia{i});
    if contains(outcome_analgesia{i}, 'rank')
        % take rank_Yint in the model
        % thermal
        mdl_heat_analgesia{i} = fitlme(analgesia_table(analgesia_table.heat > 0,:), sprintf(['%s ~ stimLvl + rank_Yint + (stimLvl + rank_Yint | family_ID) + ',...
            '(stim_mz + rank_Yint_mz + int_mz -1 | participant_ID) + ',...
            '(stim_dz + rank_Yint_dz + int_dz -1 | participant_ID)'], outcome_analgesia{i}),'FitMethod','REML');
        
        [~,~,STATS_heat_analgesia{i}] = fixedEffects(mdl_heat_analgesia{i},'dfmethod','satterthwaite');
        
        % pressure
        mdl_press_analgesia{i} = fitlme(analgesia_table(analgesia_table.heat <= 0,:), sprintf(['%s ~ stimLvl + rank_Yint + (stimLvl + rank_Yint | family_ID) + ',...
            '(stim_mz + rank_Yint_mz + int_mz -1 | participant_ID) + ',...
            '(stim_dz + rank_Yint_dz + int_dz -1 | participant_ID)'], outcome_analgesia{i}),'FitMethod','REML');
        
        [~,~,STATS_press_analgesia{i}] = fixedEffects(mdl_press_analgesia{i},'dfmethod','satterthwaite');
        
   else
        % thermal
        mdl_heat_analgesia{i} = fitlme(analgesia_table(analgesia_table.heat > 0,:), sprintf(['%s ~ stimLvl + Yint + (stimLvl + Yint | family_ID) + ',...
            '(stim_mz + Yint_mz + int_mz -1 | participant_ID) + ',...
            '(stim_dz + Yint_dz + int_dz -1 | participant_ID)'], outcome_analgesia{i}),'FitMethod','REML');
        
        [~,~,STATS_heat_analgesia{i}] = fixedEffects(mdl_heat_analgesia{i},'dfmethod','satterthwaite');
        
        % pressure
        mdl_press_analgesia{i} = fitlme(analgesia_table(analgesia_table.heat <= 0,:), sprintf(['%s ~ stimLvl + Yint + (stimLvl + Yint | family_ID) + ',...
            '(stim_mz + Yint_mz + int_mz -1 | participant_ID) + ',...
            '(stim_dz + Yint_dz + int_dz -1 | participant_ID)'], outcome_analgesia{i}),'FitMethod','REML');
        
        [~,~,STATS_press_analgesia{i}] = fixedEffects(mdl_press_analgesia{i},'dfmethod','satterthwaite');        
    end
end

%% create a table with stats
stats_table_analgesia = table;
stats_table_analgesia.outcome = repmat(outcome_analgesia',2,1);
stats_table_analgesia.modality = repelem({'thermal';'mechanical'},length(outcome_analgesia));
for ind = 1:length(STATS_heat_analgesia)
    stats_table_analgesia.estimate_yint(ind) = STATS_heat_analgesia{ind}.Estimate(3);
    stats_table_analgesia.estimate_yint(ind + length(outcome_analgesia)) = STATS_press_analgesia{ind}.Estimate(3);
    stats_table_analgesia.SE_yint(ind) = STATS_heat_analgesia{ind}.SE(3);
    stats_table_analgesia.SE_yint(ind + length(outcome_analgesia)) = STATS_press_analgesia{ind}.SE(3);
    stats_table_analgesia.t_yint(ind) = STATS_heat_analgesia{ind}.tStat(3);
    stats_table_analgesia.t_yint(ind + length(outcome_analgesia)) = STATS_press_analgesia{ind}.tStat(3);
    stats_table_analgesia.DF_yint(ind) = STATS_heat_analgesia{ind}.DF(3);
    stats_table_analgesia.DF_yint(ind + length(outcome_analgesia)) = STATS_press_analgesia{ind}.DF(3);
    stats_table_analgesia.p_yint(ind) = STATS_heat_analgesia{ind}.pValue(3);
    stats_table_analgesia.p_yint(ind + length(outcome_analgesia)) = STATS_press_analgesia{ind}.pValue(3);
    stats_table_analgesia.lower_yint(ind) = STATS_heat_analgesia{ind}.Lower(3);
    stats_table_analgesia.lower_yint(ind + length(outcome_analgesia)) = STATS_press_analgesia{ind}.Lower(3);
    stats_table_analgesia.upper_yint(ind) = STATS_heat_analgesia{ind}.Upper(3);
    stats_table_analgesia.upper_yint(ind + length(outcome_analgesia)) = STATS_press_analgesia{ind}.Upper(3);    
end
stats_table_analgesia.direction(stats_table_analgesia.estimate_yint > 0) = {'+'};
stats_table_analgesia.direction(stats_table_analgesia.estimate_yint < 0) = {'-'};
stats_table_analgesia.significance(stats_table_analgesia.p_yint >= 0.05) = {'n.s.'};
stats_table_analgesia.significance(stats_table_analgesia.p_yint < 0.05) = {'*'};
stats_table_analgesia.significance(stats_table_analgesia.p_yint < 0.01) = {'**'};
stats_table_analgesia.significance(stats_table_analgesia.p_yint < 0.001) = {'***'};

%% save model and table
if z_score
    writetable(stats_table_analgesia, fullfile(main_dir, 'stats_analgesia_models_z.csv'));
    save(fullfile(main_dir, 'models_z'), 'mdl_heat_analgesia', 'mdl_press_analgesia','STATS_heat_analgesia', 'STATS_press_analgesia', 'outcome_analgesia', "-append")
else
    writetable(stats_table_analgesia, fullfile(main_dir, 'stats_analgesia_models.csv'));
    save(fullfile(main_dir, 'models'), 'mdl_heat_analgesia', 'mdl_press_analgesia','STATS_heat_analgesia', 'STATS_press_analgesia', 'outcome_analgesia', "-append")
end

%% organize table for paper
table_paper_corr_analgesia = table;
table_paper_corr_analgesia.outcome = repelem(outcome_analgesia',2);
table_paper_corr_analgesia.modality = repmat({'thermal';'mechanical'},length(outcome_analgesia),1);
stats_table_analgesia_matched = join(table_paper_corr_analgesia, stats_table_analgesia,'Keys',{'outcome', 'modality'});
table_paper_corr_analgesia.estimate = round(stats_table_analgesia_matched.estimate_yint,3);
table_paper_corr_analgesia.CI = strcat('[', num2str(round(stats_table_analgesia_matched.lower_yint,3)), ', ', num2str(round(stats_table_analgesia_matched.upper_yint,3)), ']');
table_paper_corr_analgesia.SE = round(stats_table_analgesia_matched.SE_yint,3);
table_paper_corr_analgesia.t = round(stats_table_analgesia_matched.t_yint, 3);
table_paper_corr_analgesia.DF = round(stats_table_analgesia_matched.DF_yint,1);
table_paper_corr_analgesia.p = round(stats_table_analgesia_matched.p_yint,3);
if z_score
    writetable(table_paper_corr_analgesia, fullfile(main_dir, '/table_paper_corr_analgesia_z.csv'));
else
    writetable(table_paper_corr_analgesia, fullfile(main_dir, '/table_paper_corr_analgesia.csv'));
end