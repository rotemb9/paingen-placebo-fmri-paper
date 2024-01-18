%%% this script runs the additional models (with demographic covariates,
%%% correlations with expectations, correlations between modalities, effect
%%% sizes of NPS etc.)
%%% written by Rotem Botvinik-Nezer

%% CHANGE PATHS TO YOUR LOCAL PATHS
% add canlab tools to path
addpath(genpath('change to the local path of your copy of canlab tools'));

%%% get data
data_path = pwd; % change as needed
z_score = true; % use z scored measures
if z_score
    tt = readtable(fullfile(data_path, 'pain_evoked_brain_measures_z.csv'));
    analgesia_table = readtable(fullfile(data_path, 'pain_evoked_analgesia_z.csv'));
    tt_raw = readtable(fullfile(data_path, 'pain_evoked_brain_measures.csv'));
    analgesia_table_raw = readtable(fullfile(data_path, 'pain_evoked_analgesia.csv'));
else
    tt = readtable(fullfile(data_path, 'pain_evoked_brain_measures.csv'));
    analgesia_table = readtable(fullfile(data_path, 'pain_evoked_analgesia.csv'));
end

%%% test covariates

%% sex or age differences in placebo effects?
analgesia_table.sex = strcmp(analgesia_table.sex, 'M'); % female-0, male-1
analgesia_table.sex = analgesia_table.sex - 0.5;
analgesia_table.age_z(analgesia_table.heat > 0) = zscore(analgesia_table.age(analgesia_table.heat > 0));
analgesia_table.age_z(analgesia_table.heat <= 0) = zscore(analgesia_table.age(analgesia_table.heat <= 0));

placebo_by_demographic_thermal = fitlme(analgesia_table(analgesia_table.heat > 0,:), ['Yint ~ stimLvl + sex + age_z + (stimLvl + sex | family_ID) + ',...
            '(stim_mz + int_mz -1 | participant_ID) + ',...
            '(stim_dz + int_dz -1 | participant_ID)'],'FitMethod','REML');
        
placebo_by_demographic_mechanical = fitlme(analgesia_table(analgesia_table.heat <= 0,:), ['Yint ~ stimLvl + sex + age_z + (stimLvl + sex | family_ID) + ',...
            '(stim_mz + int_mz -1 | participant_ID) + ',...
            '(stim_dz + int_dz -1 | participant_ID)'],'FitMethod','REML');       
        
[~,~,STATS_placebo_by_demographic_thermal] = fixedEffects(placebo_by_demographic_thermal,'dfmethod','satterthwaite');
[~,~,STATS_placebo_by_demographic_mechanical] = fixedEffects(placebo_by_demographic_mechanical,'dfmethod','satterthwaite');

%%% robustness of results for pain ratings, NPS and SIIPS: controlling for
%%% sex and age
tt.sex = strcmp(tt.sex, 'M'); % female-0, male-1
tt.sex = tt.sex - 0.5;
tt.age_z(tt.heat > 0) = zscore(tt.age(tt.heat > 0));
tt.age_z(tt.heat <= 0) = zscore(tt.age(tt.heat <= 0));

% pain ratings
mdl_heat_demographic_Yint = fitlme(tt(tt.heat > 0,:), 'Yint ~ stimLvl + prodicaine + sex + age_z + (stimLvl + prodicaine + sex | family_ID) + (stim_mz + prodicaine_mz + int_mz -1 | participant_ID) + (stim_dz + prodicaine_dz + int_dz -1 | participant_ID)','FitMethod','REML');
[~,~,STATS_heat_demographic_Yint] = fixedEffects(mdl_heat_demographic_Yint,'dfmethod','satterthwaite');

mdl_press_demographic_Yint = fitlme(tt(tt.heat <= 0,:), 'Yint ~ stimLvl + prodicaine + sex + age_z + (stimLvl + prodicaine + sex | family_ID) + (stim_mz + prodicaine_mz + int_mz -1 | participant_ID) + (stim_dz + prodicaine_dz + int_dz -1 | participant_ID)','FitMethod','REML');
[~,~,STATS_press_demographic_Yint] = fixedEffects(mdl_press_demographic_Yint,'dfmethod','satterthwaite');

% NPS
mdl_heat_demographic_NPS = fitlme(tt(tt.heat > 0,:), 'nps ~ stimLvl + prodicaine + sex + age_z + (stimLvl + prodicaine + sex | family_ID) + (stim_mz + prodicaine_mz + int_mz -1 | participant_ID) + (stim_dz + prodicaine_dz + int_dz -1 | participant_ID)','FitMethod','REML');
[~,~,STATS_heat_demographic_NPS] = fixedEffects(mdl_heat_demographic_NPS,'dfmethod','satterthwaite');

mdl_press_demographic_NPS = fitlme(tt(tt.heat <= 0,:), 'nps ~ stimLvl + prodicaine + sex + age_z + (stimLvl + prodicaine + sex | family_ID) + (stim_mz + prodicaine_mz + int_mz -1 | participant_ID) + (stim_dz + prodicaine_dz + int_dz -1 | participant_ID)','FitMethod','REML');
[~,~,STATS_press_demographic_NPS] = fixedEffects(mdl_press_demographic_NPS,'dfmethod','satterthwaite');

% NPS - analgesia
analgesia_table.zygosity(:) = 1;
analgesia_table.zygosity(contains(analgesia_table.twin_ID, 'DZ') | contains(analgesia_table.twin_ID, 'OS')) = 2;
analgesia_table.Yint_mz(analgesia_table.zygosity == 1) = analgesia_table.Yint(analgesia_table.zygosity == 1);
analgesia_table.Yint_dz(analgesia_table.zygosity == 2) = analgesia_table.Yint(analgesia_table.zygosity == 2);

mdl_heat_analgesia_demographic_NPS = fitlme(analgesia_table(analgesia_table.heat > 0,:), 'nps ~ stimLvl + Yint + age_z + sex + (stimLvl + Yint + sex | family_ID) + (stim_mz + Yint_mz + int_mz -1 | participant_ID) + (stim_dz+ Yint_dz + int_dz -1 | participant_ID)','FitMethod','REML');       
[~,~,STATS_heat_analgesia_demographic_NPS] = fixedEffects(mdl_heat_analgesia_demographic_NPS,'dfmethod','satterthwaite');

mdl_press_analgesia_demographic_NPS = fitlme(analgesia_table(analgesia_table.heat <= 0,:), 'nps ~ stimLvl + Yint + age_z + sex + (stimLvl + Yint + sex | family_ID) + (stim_mz + Yint_mz + int_mz -1 | participant_ID) + (stim_dz + Yint_dz + int_dz -1 | participant_ID)','FitMethod','REML');       
[~,~,STATS_press_analgesia_demographic_NPS] = fixedEffects(mdl_press_analgesia_demographic_NPS,'dfmethod','satterthwaite');

% SIIPS
mdl_heat_demographic_SIIPS = fitlme(tt(tt.heat > 0,:), 'siips ~ stimLvl + prodicaine + sex + age_z + (stimLvl + prodicaine + sex | family_ID) + (stim_mz + prodicaine_mz + int_mz -1 | participant_ID) + (stim_dz + prodicaine_dz + int_dz -1 | participant_ID)','FitMethod','REML');
[~,~,STATS_heat_demographic_SIIPS] = fixedEffects(mdl_heat_demographic_SIIPS,'dfmethod','satterthwaite');

mdl_press_demographic_SIIPS = fitlme(tt(tt.heat <= 0,:), 'siips ~ stimLvl + prodicaine + sex + age_z + (stimLvl + prodicaine + sex | family_ID) + (stim_mz + prodicaine_mz + int_mz -1 | participant_ID) + (stim_dz + prodicaine_dz + int_dz -1 | participant_ID)','FitMethod','REML');
[~,~,STATS_press_demographic_SIIPS] = fixedEffects(mdl_press_demographic_SIIPS,'dfmethod','satterthwaite');

% SIIPS - analgesia
mdl_heat_analgesia_demographic_SIIPS = fitlme(analgesia_table(analgesia_table.heat > 0,:), 'siips ~ stimLvl + Yint + age_z + sex + (stimLvl + Yint + sex | family_ID) + (stim_mz + Yint_mz + int_mz -1 | participant_ID) + (stim_dz + Yint_dz + int_dz -1 | participant_ID)','FitMethod','REML');   
[~,~,STATS_heat_analgesia_demographic_SIIPS] = fixedEffects(mdl_heat_analgesia_demographic_SIIPS,'dfmethod','satterthwaite');

mdl_press_analgesia_demographic_SIIPS = fitlme(analgesia_table(analgesia_table.heat <= 0,:), 'siips ~ stimLvl + Yint + age_z + sex + (stimLvl + Yint + sex | family_ID) + (stim_mz + Yint_mz + int_mz -1 | participant_ID) + (stim_dz + Yint_dz + int_dz -1 | participant_ID)','FitMethod','REML');
[~,~,STATS_press_analgesia_demographic_SIIPS] = fixedEffects(mdl_press_analgesia_demographic_SIIPS,'dfmethod','satterthwaite');

%% test correlations between expectations and placebo effects
% note that values in the analgesia_table are control - placebo (so
% positive when there is analgesia)
% behavioral
expectation_data = readtable(fullfile(data_path, 'expectation_data.csv'));
fprintf('There are %d participants with valid expectation data, to be included in the models with expectation\n', height(expectation_data));
analgesia_table = outerjoin(analgesia_table, expectation_data);
analgesia_table = removevars(analgesia_table, 'participant_ID_expectation_data');
analgesia_table = renamevars(analgesia_table, {'participant_ID_analgesia_table'}, {'participant_ID'});

tt = outerjoin(tt, expectation_data);
tt = removevars(tt, 'participant_ID_expectation_data');
tt = renamevars(tt, {'participant_ID_tt', 'prodicaine_expect_post'}, {'participant_ID', 'expect'});
analgesia_table.expect_z(:) = nan;
analgesia_table.expect_z(analgesia_table.heat > 0 & ~isnan(analgesia_table.expect)) = zscore(analgesia_table.expect(analgesia_table.heat > 0 & ~isnan(analgesia_table.expect)));
analgesia_table.expect_z(analgesia_table.heat <= 0 & ~isnan(analgesia_table.expect)) = zscore(analgesia_table.expect(analgesia_table.heat <= 0 & ~isnan(analgesia_table.expect)));
tt.expect_z(:) = nan;
tt.expect_z(tt.heat > 0 & ~isnan(tt.expect)) = zscore(tt.expect(tt.heat > 0 & ~isnan(tt.expect)));
tt.expect_z(tt.heat <= 0 & ~isnan(tt.expect)) = zscore(tt.expect(tt.heat <= 0 & ~isnan(tt.expect)));

analgesia_table.expect_z_mz(analgesia_table.zygosity == 1) = analgesia_table.expect_z(analgesia_table.zygosity == 1);
analgesia_table.expect_z_dz(analgesia_table.zygosity == 2) = analgesia_table.expect_z(analgesia_table.zygosity == 2);
tt.expect_z_mz(tt.zygosity == 1) = tt.expect_z(tt.zygosity == 1);
tt.expect_z_dz(tt.zygosity == 2) = tt.expect_z(tt.zygosity == 2);

placebo_exp_heat = fitlme(analgesia_table(analgesia_table.heat > 0, :), 'Yint ~ expect_z + stimLvl + (expect_z + stimLvl | family_ID) + (stim_mz + expect_z_mz + int_mz -1 | participant_ID) + (stim_dz + expect_z_dz + int_dz -1 | participant_ID)','FitMethod','REML');
[~,~,STATS_placebo_exp_heat] = fixedEffects(placebo_exp_heat,'dfmethod','satterthwaite');
placebo_exp_press = fitlme(analgesia_table(analgesia_table.heat < 0, :), 'Yint ~ expect_z + stimLvl + (expect_z + stimLvl | family_ID) + (stim_mz + expect_z_mz + int_mz -1 | participant_ID) + (stim_dz + expect_z_dz + int_dz -1 | participant_ID)','FitMethod','REML');
[~,~,STATS_placebo_exp_press] = fixedEffects(placebo_exp_press,'dfmethod','satterthwaite');
% NPS
placebo_exp_heat_nps = fitlme(analgesia_table(analgesia_table.heat > 0, :), 'nps ~ expect_z + stimLvl + (expect_z + stimLvl | family_ID) + (stim_mz + expect_z_mz + int_mz -1 | participant_ID) + (stim_dz + expect_z_dz + int_dz -1 | participant_ID)','FitMethod','REML');
[~,~,STATS_placebo_exp_heat_nps] = fixedEffects(placebo_exp_heat_nps,'dfmethod','satterthwaite');
placebo_exp_press_nps = fitlme(analgesia_table(analgesia_table.heat < 0, :), 'nps ~ expect_z + stimLvl + (expect_z + stimLvl | family_ID) + (stim_mz + expect_z_mz + int_mz -1 | participant_ID) + (stim_dz + expect_z_dz + int_dz -1 | participant_ID)','FitMethod','REML');
[~,~,STATS_placebo_exp_press_nps] = fixedEffects(placebo_exp_press_nps,'dfmethod','satterthwaite');
% SIIPS
placebo_exp_heat_siips = fitlme(analgesia_table(analgesia_table.heat > 0, :), 'siips ~ expect_z + stimLvl + (expect_z + stimLvl | family_ID) + (stim_mz + expect_z_mz + int_mz -1 | participant_ID) + (stim_dz + expect_z_dz + int_dz -1 | participant_ID)','FitMethod','REML');
[~,~,STATS_placebo_exp_heat_siips] = fixedEffects(placebo_exp_heat_siips,'dfmethod','satterthwaite');
placebo_exp_press_siips = fitlme(analgesia_table(analgesia_table.heat < 0, :), 'siips ~ expect_z + stimLvl + (expect_z + stimLvl | family_ID) + (stim_mz + expect_z_mz + int_mz -1 | participant_ID) + (stim_dz + expect_z_dz + int_dz -1 | participant_ID)','FitMethod','REML');
[~,~,STATS_placebo_exp_press_siips] = fixedEffects(placebo_exp_press_siips,'dfmethod','satterthwaite');

%% test correlations between placebo effects in the thermal and mechanical modalities
analgesia_table(isnan(analgesia_table.stimLvl), :) = [];
modalities = {'heat', 'press'};
analgesia_table.modality(:) = modalities(1 + (analgesia_table.heat < 0));

% behavioral
placebo_heat_press_data = unstack(analgesia_table(:, {'participant_ID', 'family_ID', 'modality', 'Yint', 'zygosity', 'int_mz', 'int_dz'}), 'Yint', 'modality');
placebo_heat_press_data.heat_mz(placebo_heat_press_data.zygosity == 1) = placebo_heat_press_data.heat(placebo_heat_press_data.zygosity == 1);
placebo_heat_press_data.heat_dz(placebo_heat_press_data.zygosity == 2) = placebo_heat_press_data.heat(placebo_heat_press_data.zygosity == 2);
placebo_heat_press = fitlme(placebo_heat_press_data, 'press ~ heat + (heat | family_ID) + (heat_mz + int_mz -1 | participant_ID) + (heat_dz + int_dz -1 | participant_ID)','FitMethod','REML');
[~,~,STATS_placebo_heat_press] = fixedEffects(placebo_heat_press,'dfmethod','satterthwaite');

% NPS
placebo_heat_press_data_nps = unstack(analgesia_table(:, {'participant_ID', 'family_ID', 'modality', 'nps', 'zygosity', 'int_mz', 'int_dz'}), 'nps', 'modality');
placebo_heat_press_data_nps.heat_mz(placebo_heat_press_data_nps.zygosity == 1) = placebo_heat_press_data_nps.heat(placebo_heat_press_data_nps.zygosity == 1);
placebo_heat_press_data_nps.heat_dz(placebo_heat_press_data_nps.zygosity == 2) = placebo_heat_press_data_nps.heat(placebo_heat_press_data_nps.zygosity == 2);
placebo_heat_press_nps = fitlme(placebo_heat_press_data_nps, 'press ~ heat + (heat | family_ID) + (heat_mz + int_mz -1 | participant_ID) + (heat_dz + int_dz -1 | participant_ID)','FitMethod','REML');
[~,~,STATS_placebo_heat_press_nps] = fixedEffects(placebo_heat_press_nps,'dfmethod','satterthwaite');

% SIIPS
placebo_heat_press_data_siips = unstack(analgesia_table(:, {'participant_ID', 'family_ID', 'modality', 'siips', 'zygosity', 'int_mz', 'int_dz'}), 'siips', 'modality');
placebo_heat_press_data_siips.heat_mz(placebo_heat_press_data_siips.zygosity == 1) = placebo_heat_press_data_siips.heat(placebo_heat_press_data_siips.zygosity == 1);
placebo_heat_press_data_siips.heat_dz(placebo_heat_press_data_siips.zygosity == 2) = placebo_heat_press_data_siips.heat(placebo_heat_press_data_siips.zygosity == 2);
placebo_heat_press_siips = fitlme(placebo_heat_press_data_siips, 'press ~ heat + (heat | family_ID) + (heat_mz + int_mz -1 | participant_ID) + (heat_dz + int_dz -1 | participant_ID)','FitMethod','REML');
[~,~,STATS_placebo_heat_press_siips] = fixedEffects(placebo_heat_press_siips,'dfmethod','satterthwaite');

%% test correlations between NPS and SIIPS placebo-induced reductions
analgesia_table.siips_mz(analgesia_table.zygosity == 1) = analgesia_table.siips(analgesia_table.zygosity == 1);
analgesia_table.siips_dz(analgesia_table.zygosity == 2) = analgesia_table.siips(analgesia_table.zygosity == 2);
placebo_nps_siips_heat = fitlme(analgesia_table(analgesia_table.heat > 0, :), 'nps ~ siips + stimLvl + (siips + stimLvl | family_ID) + (stim_mz + siips_mz + int_mz -1 | participant_ID) + (stim_dz + siips_dz + int_dz -1 | participant_ID)','FitMethod','REML');
[~,~,STATS_placebo_nps_siips_heat] = fixedEffects(placebo_nps_siips_heat,'dfmethod','satterthwaite');
placebo_nps_siips_press = fitlme(analgesia_table(analgesia_table.heat < 0, :), 'nps ~ siips + stimLvl + (siips + stimLvl | family_ID) + (stim_mz + siips_mz + int_mz -1 | participant_ID) + (stim_dz + siips_dz + int_dz -1 | participant_ID)','FitMethod','REML');
[~,~,STATS_placebo_nps_siips_press] = fixedEffects(placebo_nps_siips_press,'dfmethod','satterthwaite');

%% compare placebo effects in thermal vs. mechanical modalities
analgesia_table.heat_mz(analgesia_table.zygosity == 1) = analgesia_table.heat(analgesia_table.zygosity == 1);
analgesia_table.heat_dz(analgesia_table.zygosity == 2) = analgesia_table.heat(analgesia_table.zygosity == 2);

% behavioral ratings 
analgesia_thermal_mechanical_comp = fitlme(analgesia_table, 'Yint ~ stimLvl + heat + (stimLvl + heat | family_ID) + (stim_mz + heat_mz + int_mz -1 | participant_ID) + (stim_dz + heat_dz + int_dz -1 | participant_ID)','FitMethod','REML');
[~,~,STATS_analgesia_thermal_mechanical_comp] = fixedEffects(analgesia_thermal_mechanical_comp,'dfmethod','satterthwaite');

% nps
analgesia_thermal_mechanical_comp_nps = fitlme(analgesia_table, 'nps ~ stimLvl + heat + (stimLvl + heat | family_ID) + (stim_mz + heat_mz + int_mz -1 | participant_ID) + (stim_dz + heat_dz + int_dz -1 | participant_ID)','FitMethod','REML');
[~,~,STATS_analgesia_thermal_mechanical_comp_nps] = fixedEffects(analgesia_thermal_mechanical_comp_nps,'dfmethod','satterthwaite');

% siips
analgesia_thermal_mechanical_comp_siips = fitlme(analgesia_table, 'siips ~ stimLvl + heat + (stimLvl + heat | family_ID) + (stim_mz + heat_mz + int_mz -1 | participant_ID) + (stim_dz + heat_dz + int_dz -1 | participant_ID)','FitMethod','REML');
[~,~,STATS_analgesia_thermal_mechanical_comp_siips] = fixedEffects(analgesia_thermal_mechanical_comp_siips,'dfmethod','satterthwaite');

%%% compute pain ratings, NPS and SIIPS effect sizes for pain vs. rest (baseline)
% get the mean response for each participant
% compute mean, SD and Cohen's d
sids = unique(tt_raw.participant_ID);
num_subjs = length(sids);
[Yint_mean_analgesia_heat, Yint_mean_analgesia_press, NPS_mean_heat, NPS_mean_press, NPS_mean_analgesia_heat, NPS_mean_analgesia_press, SIIPS_mean_heat, SIIPS_mean_press, SIIPS_mean_analgesia_heat, SIIPS_mean_analgesia_press] = deal(zeros(num_subjs, 1));
for sub_ind = 1:num_subjs
   cur_sid = sids(sub_ind);
   
   % Yint
   Yint_mean_analgesia_heat(sub_ind) = mean(analgesia_table_raw.Yint(strcmp(analgesia_table_raw.participant_ID, cur_sid) & analgesia_table_raw.heat > 0));  
   Yint_mean_analgesia_press(sub_ind) = mean(analgesia_table_raw.Yint(strcmp(analgesia_table_raw.participant_ID, cur_sid) & analgesia_table_raw.heat < 0));  
   
   % NPS
   NPS_mean_heat(sub_ind) = mean(tt_raw.nps(strcmp(tt_raw.participant_ID, cur_sid) & tt_raw.heat > 0));  
   NPS_mean_press(sub_ind) = mean(tt_raw.nps(strcmp(tt_raw.participant_ID, cur_sid) & tt_raw.heat < 0));
   NPS_mean_analgesia_heat(sub_ind) = mean(analgesia_table_raw.nps(strcmp(analgesia_table_raw.participant_ID, cur_sid) & analgesia_table_raw.heat > 0));  
   NPS_mean_analgesia_press(sub_ind) = mean(analgesia_table_raw.nps(strcmp(analgesia_table_raw.participant_ID, cur_sid) & analgesia_table_raw.heat < 0));  
   
   % SIIPS
   SIIPS_mean_heat(sub_ind) = mean(tt_raw.siips(strcmp(tt_raw.participant_ID, cur_sid) & tt_raw.heat > 0));  
   SIIPS_mean_press(sub_ind) = mean(tt_raw.siips(strcmp(tt_raw.participant_ID, cur_sid) & tt_raw.heat < 0)); 
   SIIPS_mean_analgesia_heat(sub_ind) = mean(analgesia_table_raw.siips(strcmp(analgesia_table_raw.participant_ID, cur_sid) & analgesia_table_raw.heat > 0));  
   SIIPS_mean_analgesia_press(sub_ind) = mean(analgesia_table_raw.siips(strcmp(analgesia_table_raw.participant_ID, cur_sid) & analgesia_table_raw.heat < 0));  
end

% Yint
yint_mean_analgesia_heat = mean(Yint_mean_analgesia_heat(~isnan(Yint_mean_analgesia_heat)));
yint_sd_analgesia_heat = std(Yint_mean_analgesia_heat(~isnan(Yint_mean_analgesia_heat)));
yint_analgesia_d_heat = yint_mean_analgesia_heat/yint_sd_analgesia_heat;
yint_mean_analgesia_press = mean(Yint_mean_analgesia_press(~isnan(Yint_mean_analgesia_press)));
yint_sd_analgesia_press = std(Yint_mean_analgesia_press(~isnan(Yint_mean_analgesia_press)));
yint_analgesia_d_press = yint_mean_analgesia_press/yint_sd_analgesia_press;

% nps
% pain vs. rest
nps_mean_heat = mean(NPS_mean_heat);
nps_sd_heat = std(NPS_mean_heat);
nps_pain_rest_d_heat = nps_mean_heat/nps_sd_heat;
nps_mean_press = mean(NPS_mean_press);
nps_sd_press = std(NPS_mean_press);
nps_pain_rest_d_press = nps_mean_press/nps_sd_press;
% placebo vs. control (analgesia)
nps_mean_analgesia_heat = mean(NPS_mean_analgesia_heat);
nps_sd_analgesia_heat = std(NPS_mean_analgesia_heat);
nps_analgesia_d_heat = nps_mean_analgesia_heat/nps_sd_analgesia_heat;
nps_mean_analgesia_press = mean(NPS_mean_analgesia_press);
nps_sd_analgesia_press = std(NPS_mean_analgesia_press);
nps_analgesia_d_press = nps_mean_analgesia_press/nps_sd_analgesia_press;

% SIIPS
% pain vs. rest
siips_mean_heat = mean(SIIPS_mean_heat);
siips_sd_heat = std(SIIPS_mean_heat);
siips_pain_rest_d_heat = siips_mean_heat/siips_sd_heat;
siips_mean_press = mean(SIIPS_mean_press);
siips_sd_press = std(SIIPS_mean_press);
siips_pain_rest_d_press = siips_mean_press/siips_sd_press;
% placebo vs. control (analgesia)
siips_mean_analgesia_heat = mean(SIIPS_mean_analgesia_heat);
siips_sd_analgesia_heat = std(SIIPS_mean_analgesia_heat);
siips_analgesia_d_heat = siips_mean_analgesia_heat/siips_sd_analgesia_heat;
siips_mean_analgesia_press = mean(SIIPS_mean_analgesia_press);
siips_sd_analgesia_press = std(SIIPS_mean_analgesia_press);
siips_analgesia_d_press = siips_mean_analgesia_press/siips_sd_analgesia_press;

