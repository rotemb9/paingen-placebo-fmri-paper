# paingen-placebo-fmri-paper
Code and (behavioral) data for the paper on placebo effects based on the Paingen dataset.
The MRI dataset is available on OpenNeuro, ds004746 (doi:10.18112/openneuro.ds004746.v1.0.0).

Since the sample includes twins, participant ID, twin ID and family ID are available in most data files. The first two letters of the twin_ID indicate the twin type (MZ = monozygotic, DZ = dizygotic, OS = dizygotic opposite sex).

Data in this repository include:
* merged_firstlvl_stim_betas_to_share.csv: mean behavioral pain ratings and response times (RTs) for each condition of each participant. Yint = intensity rating; Yunp = unpleasantness rating; stimLvl = the stimulus level, 3- low, 4- med, 5- high; heat = the modality, 1-thermal, 0-mechanical; procaine = the condition, 1=placebo, 0=control.
The corresponding nii.gz file with merged first level beta maps for each condition of each participant is available with the dataset on OpenNeuro, under derivatives, and is needed for running the script 'evaluate_pain_evoked_measures.m'.
* expectation_data: the expected Prodicaine efficacy reported by each participant (on a scale 0-100) after the conditioning phase and before the MRI scan.
* pain_evoked_brain_measures.csv: Created by the script 'evaluate_pain_evoked_measures.m', and includes the behavioral and brain measures for each condition of each participant/
* pain_evoked_brain_measures_z.csv: Same as pain_evoked_brain_measures.csv, but with values z scored within modality.
* pain_evoked_analgesia.csv: Created by the script 'evaluate_pain_evoked_measures.m', and includes the behavioral and brain measures for control minus placebo, for each participant, modality and stimulus level.
* pain_evoked_analgesia_z.csv: Same as pain_evoked_analgesia.csv, but with values z scored within modality.
* stats_condition_models: The stats from the mixed effects models, placebo, intensity and placebo X intensity effects.
* stats_condition_models_z: Same as stats_condition_models, but with models that are based on z scored values (and thus with normalized estimates).
* stats_analgesia_models: The stats from the mixed effects models, correlations between behavioral analgesia and neural placebo-induced reductions.
* stats_analgesia_models_z: Same as stats_analgesia_models, but with models that are based on z scored values (and thus with normalized estimates).
* table_paper_*.csv: the stats organized in a table as in the paper.

Scripts in this repository include:
* evaluate_pain_evoked_measures.m: Uses CANlab imaging analysis tools (https://canlab.github.io/) to compute the score of each brain signature and ROI for each condition of each participant.
* pain_evoked_models.m: Uses the measures estimated with evaluate_pain_evoked_measures.m (and saved in pain_evoked_brain_measures_z.csv) to estimate the placebo, intensity, and placebo X intensity effects for the behavioral and the neural scores.
* pain_evoked_models_ind_diff.m: Uses the measures estimated with evaluate_pain_evoked_measures.m (and saved in pain_evoked_analgesia_z.csv) to estimate the correlations between behavioral analgesia and neural placebo-induced reductions.
* additional_models.m: Uses the measures estimated with evaluate_pain_evoked_measures.m and runs additional models, including covariates, correlations between modalities, differences between modalities, effect sizes, and more.
* BF_analysis.R: Runs and saves the results of the Bayes Factor analyses.
* plots.R: Creates all the figures of the paper. The last part also computes some of the condition means that are reported in the paper.

