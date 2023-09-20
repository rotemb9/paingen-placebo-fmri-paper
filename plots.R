# This script creates plots for the paingen fmri placebo paper
# Written by Rotem Botvinik-Nezer

################
# PREPERATIONS #
################

# load libraries
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(Rmisc)
library(Hmisc)
library(png, include.only = 'readPNG')
library(gridExtra)

## Clear workspace
rm(list=ls())

### SETTINGS
fig_output_path = "figures"

### read data
data = read.csv("pain_evoked_brain_measures.csv")
data_z = read.csv("pain_evoked_brain_measures_z.csv")
data_analgesia = read.csv("pain_evoked_analgesia.csv")
data_analgesia_z = read.csv("pain_evoked_analgesia_z.csv")
stats_condition_z = read.csv("table_paper_condition_z.csv")
stats_analgesia_z = read.csv("table_paper_corr_analgesia_z.csv")

## prepare labels for plots
stimLvl_labels = c("low", "med", "high")
prodicaine_labels = c("control", "placebo")
heat_labels = c("mechanical", "thermal")

data$stimLvl_labels = stimLvl_labels[data$stimLvl * 2 + 2]
data$stimLvl_labels = factor(data$stimLvl_labels, ordered = TRUE, levels = c("low", "med", "high"))
data$prodicaine_labels = prodicaine_labels[data$prodicaine + 1.5]
data$heat_labels = heat_labels[data$heat + 1.5]
data$heat_labels = factor(data$heat_labels, levels = c("thermal", "mechanical"))
data_z$stimLvl_labels = stimLvl_labels[data_z$stimLvl * 2 + 2]
data_z$stimLvl_labels = factor(data_z$stimLvl_labels, ordered = TRUE, levels = c("low", "med", "high"))
data_z$prodicaine_labels = prodicaine_labels[data_z$prodicaine + 1.5]
data_z$heat_labels = heat_labels[data_z$heat + 1.5]
data_z$heat_labels = factor(data_z$heat_labels, levels = c("thermal", "mechanical"))

data_analgesia$stimLvl_labels = stimLvl_labels[data_analgesia$stimLvl * 2 + 2]
data_analgesia$stimLvl_labels = factor(data_analgesia$stimLvl_labels, ordered = TRUE, levels = c("low", "med", "high"))
data_analgesia$heat_labels = heat_labels[data_analgesia$heat + 1.5]
data_analgesia$heat_labels = factor(data_analgesia$heat_labels, levels = c("thermal", "mechanical"))
data_analgesia_z$stimLvl_labels = stimLvl_labels[data_analgesia_z$stimLvl * 2 + 2]
data_analgesia_z$stimLvl_labels = factor(data_analgesia_z$stimLvl_labels, ordered = TRUE, levels = c("low", "med", "high"))
data_analgesia_z$heat_labels = heat_labels[data_analgesia_z$heat + 1.5]
data_analgesia_z$heat_labels = factor(data_analgesia_z$heat_labels, levels = c("thermal", "mechanical"))

## convert categorical variables into factors
data$participant_ID = as.factor(data$participant_ID)
data$family_ID = as.factor(data$family_ID)
data$stimLvl = as.factor(data$stimLvl)
data$prodicaine = as.factor(data$prodicaine)
data$heat = as.factor(data$heat)

data_z$participant_ID = as.factor(data_z$participant_ID)
data_z$family_ID = as.factor(data_z$family_ID)
data_z$stimLvl = as.factor(data_z$stimLvl)
data_z$prodicaine = as.factor(data_z$prodicaine)
data_z$heat = as.factor(data_z$heat)

data_analgesia$participant_ID = as.factor(data_analgesia$participant_ID)
data_analgesia$family_ID = as.factor(data_analgesia$family_ID)
data_analgesia$stimLvl = as.factor(data_analgesia$stimLvl)
data_analgesia$heat = as.factor(data_analgesia$heat)

data_analgesia_z$participant_ID = as.factor(data_analgesia_z$participant_ID)
data_analgesia_z$family_ID = as.factor(data_analgesia_z$family_ID)
data_analgesia_z$stimLvl = as.factor(data_analgesia_z$stimLvl)
data_analgesia_z$heat = as.factor(data_analgesia_z$heat)

# prepare asterisks for the significance level of placebo and stim level effects
stats_condition_z$asterisks = ""
stats_condition_z$asterisks[stats_condition_z$p < 0.05] = "*"
stats_condition_z$asterisks[stats_condition_z$p < 0.01] = "**"
stats_condition_z$asterisks[stats_condition_z$p < 0.001] = "***"
stats_condition_z$effect[stats_condition_z$effect == "Stim level"] = "stim_level"
stats_condition_z$effect[stats_condition_z$effect == "Placebo"] = "placebo"
stats_condition_z = dplyr::rename(stats_condition_z, measure = outcome)

stats_analgesia_z$asterisks = ""
stats_analgesia_z$asterisks[stats_analgesia_z$p < 0.05] = "*"
stats_analgesia_z$asterisks[stats_analgesia_z$p < 0.01] = "**"
stats_analgesia_z$asterisks[stats_analgesia_z$p < 0.001] = "***"
stats_analgesia_z = dplyr::rename(stats_analgesia_z, measure = outcome)

#########
# PLOTS #
#########
## create summarized df for plotting the behavioral ratings and brain scores (first NPS and SIIPS)
data_long = gather(select(data, c(participant_ID, family_ID, stimLvl_labels, prodicaine_labels, heat_labels, Yint, nps, siips)), measure, score, c(Yint, nps, siips), factor_key=TRUE)
data_for_plot = aggregate(data_long$score,list(measure = data_long$measure, modality = data_long$heat_labels, condition = data_long$prodicaine_labels, stimulus_level = data_long$stimLvl_labels, participant_ID = data_long$participant_ID), function(x) c(mean = mean(x, na.rm=TRUE)))
data_for_plot = aggregate(data_for_plot$x,list(measure = data_for_plot$measure, modality = data_for_plot$modality, condition = data_for_plot$condition, stimulus_level = data_for_plot$stimulus_level), function(x) c(mean = mean(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE), n = sum(!is.na(x))))
data_for_plot = cbind(data_for_plot[,1:ncol(data_for_plot)-1], as.data.frame(data_for_plot$x))
data_for_plot$se = data_for_plot$sd/sqrt(data_for_plot$n)

# New facet label names for modality variable
modality.labs = c("Thermal", "Mechanical")
names(modality.labs) = c("thermal", "mechanical")
# New facet label names for measure variable
measure.labs = c("Pain ratings", "NPS", "SIIPS")
names(measure.labs) = c("Yint", "nps", "siips")

stats_condition_z_to_merge = stats_condition_z[stats_condition_z$measure %in% c("nps", "siips", "Yint") & stats_condition_z$effect == "placebo", c("measure", "modality", "asterisks")]
data_for_plot = merge(data_for_plot, stats_condition_z_to_merge)
data_for_plot$asterisks[data_for_plot$condition != "control" | data_for_plot$stimulus_level != "med"] = NA;

# compute within participant SEMs for the plots, with the method based on Morey 2008 (values adjusted for inter-subject variability)
data_for_plot_behav_adjusted = summarySEwithin(data,  measurevar = "Yint", withinvars = c("stimLvl_labels", "prodicaine_labels", "heat_labels"), idvar = "participant_ID", na.rm=TRUE, conf.interval=.95)
data_for_plot_behav_adjusted = dplyr::rename(data_for_plot_behav_adjusted, c(condition = prodicaine_labels, modality = heat_labels, stimulus_level = stimLvl_labels, mean_within = Yint, se_within = se, sd_within = sd, ci_within = ci, N_within = N))
data_for_plot_Yint = merge(data_for_plot[data_for_plot$measure %in% c("Yint"),], data_for_plot_behav_adjusted)

## plot behavioral
pain_breaks = c(0.014, 0.061, 0.172, 0.354)
pain_breaks_labels = c("barely detectable", "weak", "moderate", "strong")

lineplot_behav = ggplot(data_for_plot_Yint, aes(x=stimulus_level, y=mean, colour=condition, group = condition)) + 
  geom_hline(yintercept = pain_breaks, linetype = "dashed") +
  geom_errorbar(aes(ymin=mean-se_within, ymax=mean+se_within), width=.1) +
  geom_line() +
  geom_point() +
  geom_text(aes(y = max(mean)+0.02, label=asterisks), size = 6, color = "black") +
  labs(x = "Stimulus level", y = "Mean pain rating", colour = "Condition") +
  scale_y_continuous(limits = c(0, 0.4), breaks = seq(0, 0.4, 0.05), sec.axis = sec_axis(~ ., breaks = pain_breaks, labels = pain_breaks_labels)) +
  facet_grid(. ~ modality, scales = "free", labeller = labeller(modality = modality.labs)) +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="bottom")

## plot NPS
data_for_plot_nps_adjusted = summarySEwithin(data,  measurevar = "nps", withinvars = c("stimLvl_labels", "prodicaine_labels", "heat_labels"), idvar = "participant_ID", na.rm=TRUE, conf.interval=.95)
data_for_plot_nps_adjusted = dplyr::rename(data_for_plot_nps_adjusted, c(condition = prodicaine_labels, modality = heat_labels, stimulus_level = stimLvl_labels, mean_within = nps, se_within = se, sd_within = sd, ci_within = ci, N_within = N))
data_for_plot_nps = merge(data_for_plot[data_for_plot$measure %in% c("nps"),], data_for_plot_nps_adjusted)

lineplot_nps_ver = ggplot(data_for_plot_nps, aes(x=stimulus_level, y=mean, colour=condition, group = condition)) + 
  geom_errorbar(aes(ymin=mean-se_within, ymax=mean+se_within), width=.1) +
  geom_line() +
  geom_point() +
  geom_text(aes(y = max(mean)+0.2, label=asterisks), size = 6, color = "black") +
  labs(x = "Stimulus level", y = "Mean NPS score", colour = "Condition") +
  facet_grid(modality ~ ., scales = "free", labeller = labeller(modality = modality.labs)) +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10), legend.position="bottom")

## plot SIIPS
data_for_plot_siips_adjusted = summarySEwithin(data,  measurevar = "siips", withinvars = c("stimLvl_labels", "prodicaine_labels", "heat_labels"), idvar = "participant_ID", na.rm=TRUE, conf.interval=.95)
data_for_plot_siips_adjusted = dplyr::rename(data_for_plot_siips_adjusted, c(condition = prodicaine_labels, modality = heat_labels, stimulus_level = stimLvl_labels, mean_within = siips, se_within = se, sd_within = sd, ci_within = ci, N_within = N))
data_for_plot_siips = merge(data_for_plot[data_for_plot$measure %in% c("siips"),], data_for_plot_siips_adjusted)

lineplot_siips_hor = ggplot(data_for_plot_siips, aes(x=stimulus_level, y=mean, colour=condition, group = condition)) + 
  geom_errorbar(aes(ymin=mean-se_within, ymax=mean+se_within), width=.1) +
  geom_line() +
  geom_point() +
  geom_text(aes(y = max(mean)+2, label=asterisks), size = 6, color = "black") +
  labs(x = "Stimulus level", y = "Mean SIIPS score", colour = "Condition") +
  facet_grid(. ~ modality, scales = "free", labeller = labeller(modality = modality.labs)) +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10))

lineplot_siips_ver = ggplot(data_for_plot_siips, aes(x=stimulus_level, y=mean, colour=condition, group = condition)) + 
  geom_errorbar(aes(ymin=mean-se_within, ymax=mean+se_within), width=.1) +
  geom_line() +
  geom_point() +
  geom_text(aes(y = max(mean)+2, label=asterisks), size = 6, color = "black") +
  labs(x = "Stimulus level", y = "Mean SIIPS score", colour = "Condition") +
  facet_grid(modality ~ ., scales = "free", labeller = labeller(modality = modality.labs)) +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10), legend.position = "bottom")

## bar plots with all dots
data_for_points = data_long
data_for_points = dplyr::rename(data_for_points, c(condition = prodicaine_labels, modality = heat_labels, stimulus_level = stimLvl_labels))

## behavioral
bar_dots_behav_effect_plot = ggplot(data_for_plot_Yint, aes(x=stimulus_level, y=mean, fill = condition)) + 
  geom_bar(width=0.8,colour="black",position=position_dodge(0.8), stat="identity", alpha=0.1) + # Bar plot
  geom_errorbar(position=position_dodge(0.8), width=1/4, aes(ymin=mean-se_within, ymax=mean+se_within)) +
  geom_point(data=data_for_points[data_for_points$measure %in% c("Yint"),], aes(y = score, x=stimulus_level, color=condition), position=position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8), alpha=0.1, size = 0.8) +
  geom_text(aes(y = 0.8, label=asterisks), size = 6, color = "black") +
  labs(x = "Stimulus level", y = "Pain rating", fill = "Condition") +
  facet_grid(~ modality, scales = "free", labeller = labeller(modality = modality.labs)) +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10)) +
  guides(color = "none")

## NPS + a-priori nociceptive regions
data_nociception_long = gather(select(data, c(participant_ID, family_ID, stimLvl_labels, prodicaine_labels, heat_labels, nps, starts_with("pain_pathways") | starts_with("canlab2018_mean_Bstem_PAG"))), measure, score, c(nps, starts_with("pain_pathways") | starts_with("canlab2018_mean_Bstem_PAG")), factor_key=FALSE)
nociceptive_regions = c("pain_pathways_Thal_VPLM_R", "pain_pathways_Thal_VPLM_L", "pain_pathways_Thal_MD", "pain_pathways_dpIns_L", "pain_pathways_dpIns_R", "pain_pathways_aMCC_MPFC", "canlab2018_mean_Bstem_PAG")
nociceptive_regions.labs = c("VPL/M Thal R", "VPL/M Thal L", "Med Thal", "dpIns L", "dpIns R", "aMCC", "PAG")
names(nociceptive_regions.labs) = nociceptive_regions

# condition data for nociceptive regions
for (region_ind in 1:length(nociceptive_regions)) {
  data_nociception_long$measure[data_nociception_long$measure == nociceptive_regions[region_ind]] = nociceptive_regions.labs[region_ind]
  stats_condition_z$measure[stats_condition_z$measure == nociceptive_regions[region_ind]] = nociceptive_regions.labs[region_ind]
}
data_nociception_long$measure[data_nociception_long$measure == "nps"] = "NPS"
stats_condition_z$measure[stats_condition_z$measure == "nps"] = "NPS"
data_nociception_long$measure = factor(data_nociception_long$measure)
data_nociception_long$measure = relevel(data_nociception_long$measure, ref = "NPS")
data_nociception_for_plot = aggregate(data_nociception_long$score,list(measure = data_nociception_long$measure, modality = data_nociception_long$heat_labels, condition = data_nociception_long$prodicaine_labels, stimulus_level = data_nociception_long$stimLvl_labels, participant_ID = data_nociception_long$participant_ID), function(x) c(mean = mean(x, na.rm=TRUE)))
data_nociception_for_plot = aggregate(data_nociception_for_plot$x,list(measure = data_nociception_for_plot$measure, modality = data_nociception_for_plot$modality, condition = data_nociception_for_plot$condition, stimulus_level = data_nociception_for_plot$stimulus_level), function(x) c(mean = mean(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE), n = sum(!is.na(x))))
data_nociception_for_plot = cbind(data_nociception_for_plot[,1:ncol(data_nociception_for_plot)-1], as.data.frame(data_nociception_for_plot$x))
data_nociception_for_plot$se = data_nociception_for_plot$sd/sqrt(data_nociception_for_plot$n)
data_nociception_for_points = data_nociception_long
data_nociception_for_points = dplyr::rename(data_nociception_for_points, c(condition = prodicaine_labels, modality = heat_labels, stimulus_level = stimLvl_labels))

stats_condition_z_to_merge = stats_condition_z[stats_condition_z$measure %in% nociceptive_regions.labs & stats_condition_z$effect == "placebo", c("measure", "modality", "asterisks")]
data_nociception_for_plot = merge(data_nociception_for_plot, stats_condition_z_to_merge)
data_nociception_for_plot$asterisks[data_nociception_for_plot$condition != "control" | data_nociception_for_plot$stimulus_level != "med"] = NA;

# compute the adjusted within-subject variance for each region separately, otherwise the variance increases because of variance across regions
for (region_ind in 1:length(nociceptive_regions)) {
  cur_values = summarySEwithin(data_nociception_long[data_nociception_long$measure == nociceptive_regions.labs[region_ind],],  measurevar = "score", withinvars = c("stimLvl_labels", "prodicaine_labels", "heat_labels", "measure"), idvar = "participant_ID", na.rm=TRUE, conf.interval=.95)
  if (region_ind > 1) {
    data_for_plot_nociception_adjusted = rbind(data_for_plot_nociception_adjusted, cur_values)
  } else {
    data_for_plot_nociception_adjusted = cur_values
  }
}
data_for_plot_nociception_adjusted = dplyr::rename(data_for_plot_nociception_adjusted, c(condition = prodicaine_labels, modality = heat_labels, stimulus_level = stimLvl_labels, mean_within = score, se_within = se, sd_within = sd, ci_within = ci, N_within = N))
data_for_plot_nociception = merge(data_nociception_for_plot, data_for_plot_nociception_adjusted)

# line plot, nociceptive ROIs, no NPS
lineplot_nociceptive_regions = ggplot(data_for_plot_nociception, aes(x=stimulus_level, y=mean, colour=condition, group = condition)) + 
  geom_errorbar(aes(ymin=mean-se_within, ymax=mean+se_within), width=.1) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_text(aes(y = max(mean)+ 0.01, label=asterisks), size = 6, color = "black") +
  scale_y_continuous(breaks = seq(-0.2, 0.2, 0.1), limits = c(-0.2, 0.2)) +
  labs(x = "Stimulus level", y = "Mean score", colour = "Condition") +
  facet_grid(modality ~ measure, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10))

## bars and dots - take only 95% range
# bar plot with dots - ONLY 95% / 90% RANGE OF VALUES
data_nociception_for_points_no_nps = data_nociception_for_points[data_nociception_for_points$measure != "NPS",]
for(measure in unique(data_nociception_for_points_no_nps$measure)) {
  for (stimulus_level in unique(data_nociception_for_points_no_nps$stimulus_level)) {
    for (condition in unique(data_nociception_for_points_no_nps$condition)) {
      for (modality in unique(data_nociception_for_points_no_nps$modality)) {
        cur_data = data_nociception_for_points_no_nps$score[data_nociception_for_points_no_nps$measure == measure & data_nociception_for_points_no_nps$stimulus_level == stimulus_level & data_nociception_for_points_no_nps$condition == condition & data_nociception_for_points_no_nps$modality == modality]
        ranges = quantile(cur_data, probs=c(.025, .05, .95, .975), na.rm = TRUE)
        data_nociception_for_points_no_nps$min_range[data_nociception_for_points_no_nps$measure == measure & data_nociception_for_points_no_nps$stimulus_level == stimulus_level & data_nociception_for_points_no_nps$condition == condition & data_nociception_for_points_no_nps$modality == modality] = ranges[1]
        data_nociception_for_points_no_nps$range5per[data_nociception_for_points_no_nps$measure == measure & data_nociception_for_points_no_nps$stimulus_level == stimulus_level & data_nociception_for_points_no_nps$condition == condition & data_nociception_for_points_no_nps$modality == modality] = ranges[2]
        data_nociception_for_points_no_nps$range95per[data_nociception_for_points_no_nps$measure == measure & data_nociception_for_points_no_nps$stimulus_level == stimulus_level & data_nociception_for_points_no_nps$condition == condition & data_nociception_for_points_no_nps$modality == modality] = ranges[3]
        data_nociception_for_points_no_nps$max_range[data_nociception_for_points_no_nps$measure == measure & data_nociception_for_points_no_nps$stimulus_level == stimulus_level & data_nociception_for_points_no_nps$condition == condition & data_nociception_for_points_no_nps$modality == modality] = ranges[4]
      }
    }
  }
}
data_nociception_for_points_no_nps_range90 = data_nociception_for_points_no_nps[data_nociception_for_points_no_nps$score >= data_nociception_for_points_no_nps$range5per & data_nociception_for_points_no_nps$score <= data_nociception_for_points_no_nps$range95per,]
data_nociception_for_points_no_nps_range90_1 = data_nociception_for_points_no_nps_range90[data_nociception_for_points_no_nps_range90$measure %in% nociceptive_regions.labs[1:ceil(length(nociceptive_regions.labs)/2)],]
data_nociception_for_points_no_nps_range90_2 = data_nociception_for_points_no_nps_range90[data_nociception_for_points_no_nps_range90$measure %in% nociceptive_regions.labs[ceil(length(nociceptive_regions.labs)/2) + 1:length(nociceptive_regions.labs)],]

bar_dots_nociceptive_regions_plot_range90 = ggplot(data_for_plot_nociception, aes(x=stimulus_level, y=mean, fill = condition)) + 
  geom_point(data=data_nociception_for_points_no_nps_range90, aes(y = score, x=stimulus_level, color=condition), position=position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8), alpha=0.05, size = 0.8) +
  geom_bar(width=0.8,colour="black",position=position_dodge(0.8), stat="identity", alpha=0.1) + # Bar plot
  geom_errorbar(position=position_dodge(0.8), width=1/4, aes(ymin=mean-se_within, ymax=mean+se_within)) +
  geom_text(aes(y = max(mean)+ 0.3, label=asterisks), size = 6, color = "black") +
  labs(x = "Stimulus level", y = "Score", fill = "Condition") +
  facet_grid(measure ~ modality, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10)) +
  guides(color = "none")

bar_dots_nociceptive_regions_plot_range90_1 = ggplot(data_for_plot_nociception[data_for_plot_nociception$measure %in% nociceptive_regions.labs[1:ceil(length(nociceptive_regions.labs)/2)],], aes(x=stimulus_level, y=mean, fill = condition)) + 
  geom_point(data=data_nociception_for_points_no_nps_range90_1, aes(y = score, x=stimulus_level, color=condition), position=position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8), alpha=0.05, size = 0.8) +
  geom_bar(width=0.8,colour="black",position=position_dodge(0.8), stat="identity", alpha=0.1) + # Bar plot
  geom_errorbar(position=position_dodge(0.8), width=1/4, aes(ymin=mean-se_within, ymax=mean+se_within)) +
  geom_text(aes(y = max(mean)+ 0.3, label=asterisks), size = 6, color = "black") +
  labs(x = "Stimulus level", y = "Score", fill = "Condition") +
  facet_grid(measure ~ modality, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10)) +
  guides(color = "none")

bar_dots_nociceptive_regions_plot_range90_2 = ggplot(data_for_plot_nociception[data_for_plot_nociception$measure %in% nociceptive_regions.labs[ceil(length(nociceptive_regions.labs)/2) + 1:length(nociceptive_regions.labs)],], aes(x=stimulus_level, y=mean, fill = condition)) + 
  geom_point(data=data_nociception_for_points_no_nps_range90_2, aes(y = score, x=stimulus_level, color=condition), position=position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8), alpha=0.05, size = 0.8) +
  geom_bar(width=0.8,colour="black",position=position_dodge(0.8), stat="identity", alpha=0.1) + # Bar plot
  geom_errorbar(position=position_dodge(0.8), width=1/4, aes(ymin=mean-se_within, ymax=mean+se_within)) +
  geom_text(aes(y = max(mean)+ 0.3, label=asterisks), size = 6, color = "black") +
  labs(x = "Stimulus level", y = "Score", fill = "Condition") +
  facet_grid(measure ~ modality, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10)) +
  guides(color = "none")

bar_dots_nociceptive_regions_plot_range90_combined = ggarrange(bar_dots_nociceptive_regions_plot_range90_1, bar_dots_nociceptive_regions_plot_range90_2,
                                                               labels = c("", ""),
                                                               ncol = 2, nrow = 1, 
                                                               common.legend = TRUE, legend = "bottom")


# only NPS
NPS_plot = ggplot(data_for_plot_nps, aes(x=stimulus_level, y=mean, fill = condition)) + 
  geom_point(data=data_nociception_for_points[data_nociception_for_points$measure == "NPS",], aes(y = score, x=stimulus_level, color=condition), position=position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8), alpha=0.1, size = 0.8) +
  geom_bar(width=0.8,colour="black",position=position_dodge(0.8), stat="identity", alpha=0.1) + # Bar plot
  geom_errorbar(position=position_dodge(0.8), width=1/4, aes(ymin=mean-se_within, ymax=mean+se_within)) +
  geom_text(aes(y = max(mean) + 2*sd, label=asterisks), size = 6, color = "black") +
  labs(x = "Stimulus level", y = "NPS score", color = "Condition", fill = "Condition") +
  facet_grid(. ~ modality, scales = "free", labeller = labeller(modality = modality.labs)) +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10)) +
  guides(color = "none", fill = guide_legend(title.position="top", title.hjust = 0.5))

# only NPS - ONLY 95% RANGE OF VALUES
data_nociception_for_plot_nps = data_nociception_for_plot[data_nociception_for_plot$measure == "NPS",]
data_nociception_for_points_nps = data_nociception_for_points[data_nociception_for_points$measure == "NPS",]
for (stimulus_level in unique(data_nociception_for_points_nps$stimulus_level)) {
  for (condition in unique(data_nociception_for_points_nps$condition)) {
    for (modality in unique(data_nociception_for_points_nps$modality)) {
      cur_data = data_nociception_for_points_nps$score[data_nociception_for_points_nps$stimulus_level == stimulus_level & data_nociception_for_points_nps$condition == condition & data_nociception_for_points_nps$modality == modality]
      range95 = quantile(cur_data, probs=c(.025, .975), na.rm = TRUE)
      data_nociception_for_points_nps$min_range[data_nociception_for_points_nps$stimulus_level == stimulus_level & data_nociception_for_points_nps$condition == condition & data_nociception_for_points_nps$modality == modality] = range95[1]
      data_nociception_for_points_nps$max_range[data_nociception_for_points_nps$stimulus_level == stimulus_level & data_nociception_for_points_nps$condition == condition & data_nociception_for_points_nps$modality == modality] = range95[2]
    }
  }
}
data_nociception_for_points_nps_range95 = data_nociception_for_points_nps[data_nociception_for_points_nps$score >= data_nociception_for_points_nps$min_range & data_nociception_for_points_nps$score <= data_nociception_for_points_nps$max_range,]
NPS_plot_95 = ggplot(data_for_plot_nps, aes(x=stimulus_level, y=mean, fill = condition)) + 
  geom_point(data=data_nociception_for_points_nps_range95, aes(y = score, x=stimulus_level, color=condition), position=position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8), alpha=0.1, size = 0.8) +
  geom_bar(width=0.8,colour="black",position=position_dodge(0.8), stat="identity", alpha=0.1) + # Bar plot
  geom_errorbar(position=position_dodge(0.8), width=1/4, aes(ymin=mean-se_within, ymax=mean+se_within)) +
  geom_text(aes(y = max(mean) + 2*sd, label=asterisks), size = 6, color = "black") +
  labs(x = "Stimulus level", y = "NPS score", color = "Condition", fill = "Condition") +
  facet_grid(. ~ modality, scales = "free", labeller = labeller(modality = modality.labs)) +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10)) +
  guides(color = "none", fill = guide_legend(title.position="top", title.hjust = 0.5))

### plots for SIIPS, SIIPS's subregions and all the higher level a-priori ROIs
siips_subregions = c("siips", "siipspos_dmPFC", "siipspos_R_MTG", "siipsneg_L_NAc", "siipsneg_R_LG", "siipsneg_L_STG", "siipsneg_R_TP", "siipsneg_L_ITG", "siipsneg_mid_precen")
siips_subregions.labs = c("SIIPS", "SIIPS Pos - dmPFC", "SIIPS Pos - MTG R", "SIIPS Neg - NAc L", "SIIPS Neg - LG R", "SIIPS Neg - STG L", "SIIPS Neg - TP R", "SIIPS Neg - ITG L", "SIIPS Neg - Mid Precen")
names(siips_subregions.labs) = siips_subregions

# condition data for SIIPS subregions
data_siips_subregions_long = gather(select(data, c(participant_ID, family_ID, stimLvl_labels, prodicaine_labels, heat_labels, siips, starts_with("siipspos"), starts_with("siipsneg"))), measure, score, c(siips, starts_with("siipspos"), starts_with("siipsneg")), factor_key=FALSE)
for (region_ind in 1:length(siips_subregions)) {
  data_siips_subregions_long$measure[data_siips_subregions_long$measure == siips_subregions[region_ind]] = siips_subregions.labs[region_ind]
  stats_condition_z$measure[stats_condition_z$measure == siips_subregions[region_ind]] = siips_subregions.labs[region_ind]
}
data_siips_subregions_long = data_siips_subregions_long[data_siips_subregions_long$measure %in% siips_subregions.labs,]
data_siips_subregions_long$measure = factor(data_siips_subregions_long$measure, levels = siips_subregions.labs)
data_siips_subregions_for_plot = aggregate(data_siips_subregions_long$score,list(measure = data_siips_subregions_long$measure, modality = data_siips_subregions_long$heat_labels, condition = data_siips_subregions_long$prodicaine_labels, stimulus_level = data_siips_subregions_long$stimLvl_labels), function(x) c(mean = mean(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE), n = sum(!is.na(x))))
data_siips_subregions_for_plot = cbind(data_siips_subregions_for_plot[,1:ncol(data_siips_subregions_for_plot)-1], as.data.frame(data_siips_subregions_for_plot$x))
data_siips_subregions_for_plot$se = data_siips_subregions_for_plot$sd/sqrt(data_siips_subregions_for_plot$n)
data_siips_subregions_points = data_siips_subregions_long
data_siips_subregions_points = dplyr::rename(data_siips_subregions_points, c(condition = prodicaine_labels, modality = heat_labels, stimulus_level = stimLvl_labels))

# compute the adjusted within-subject variance for each region separately, otherwise the variance increases because of variance across regions
for (region_ind in 1:length(siips_subregions)) {
  cur_values = summarySEwithin(data_siips_subregions_long[data_siips_subregions_long$measure == siips_subregions.labs[region_ind],],  measurevar = "score", withinvars = c("stimLvl_labels", "prodicaine_labels", "heat_labels", "measure"), idvar = "participant_ID", na.rm=TRUE, conf.interval=.95)
  if (region_ind > 1) {
    data_for_plot_siips_subregions_adjusted = rbind(data_for_plot_siips_subregions_adjusted, cur_values)
  } else {
    data_for_plot_siips_subregions_adjusted = cur_values
  }
}
data_for_plot_siips_subregions_adjusted = dplyr::rename(data_for_plot_siips_subregions_adjusted, c(condition = prodicaine_labels, modality = heat_labels, stimulus_level = stimLvl_labels, mean_within = score, se_within = se, sd_within = sd, ci_within = ci, N_within = N))
data_siips_subregions_for_plot = merge(data_siips_subregions_for_plot, data_for_plot_siips_subregions_adjusted)
data_siips_subregions_for_plot_no_siips = data_siips_subregions_for_plot[data_siips_subregions_for_plot$measure != "SIIPS",]

# get the significance (asterisks) info
stats_condition_z_to_merge = stats_condition_z[stats_condition_z$measure %in% siips_subregions.labs & stats_condition_z$effect == "placebo", c("measure", "modality", "asterisks")]
stats_condition_z_to_merge = stats_condition_z_to_merge[stats_condition_z_to_merge$measure != "SIIPS",]
data_siips_subregions_for_plot_no_siips = merge(data_siips_subregions_for_plot_no_siips, stats_condition_z_to_merge)
data_siips_subregions_for_plot_no_siips$asterisks[data_siips_subregions_for_plot_no_siips$condition != "control" | data_siips_subregions_for_plot_no_siips$stimulus_level != "med"] = NA;

# divide siips subregions to 2 to facilitate better layout of the plots
data_siips_subregions_for_plot_no_siips1 = data_siips_subregions_for_plot_no_siips[data_siips_subregions_for_plot_no_siips$measure %in% siips_subregions.labs[2:5],]
data_siips_subregions_for_plot_no_siips2 = data_siips_subregions_for_plot_no_siips[data_siips_subregions_for_plot_no_siips$measure %in% siips_subregions.labs[6:9],]
data_siips_subregions_points_no_siips = data_siips_subregions_points[data_siips_subregions_points$measure != "SIIPS",]
data_siips_subregions_points_no_siips1 = data_siips_subregions_points_no_siips[data_siips_subregions_points_no_siips$measure %in% siips_subregions.labs[2:5],]
data_siips_subregions_points_no_siips2 = data_siips_subregions_points_no_siips[data_siips_subregions_points_no_siips$measure %in% siips_subregions.labs[6:9],]

# SIIPS, bar with dots
SIIPS_plot = ggplot(data_for_plot_siips, aes(x=stimulus_level, y=mean, fill = condition)) + 
  geom_point(data=data_siips_subregions_points[data_siips_subregions_points$measure == "SIIPS",], aes(y = score, x=stimulus_level, color=condition), position=position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8), alpha=0.1, size = 0.8) +
  geom_bar(width=0.8,colour="black",position=position_dodge(0.8), stat="identity", alpha=0.1) + # Bar plot
  geom_errorbar(position=position_dodge(0.8), width=1/4, aes(ymin=mean-se_within, ymax=mean+se_within)) +
  geom_text(aes(y = max(mean) + 3 * sd, label=asterisks), size = 6, color = "black") +
  labs(x = "Stimulus level", y = "SIIPS score", color = "Condition", fill = "Condition") +
  facet_grid(. ~ modality, scales = "free", labeller = labeller(modality = modality.labs)) +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10)) +
  guides(color = "none", fill = guide_legend(title.position="top", title.hjust = 0.5))

# SIIPS - ONLY 95% RANGE OF VALUES
data_siips_subregions_points_siips = data_siips_subregions_points[data_siips_subregions_points$measure == "SIIPS",]
for (stimulus_level in unique(data_siips_subregions_points_siips$stimulus_level)) {
  for (condition in unique(data_siips_subregions_points_siips$condition)) {
    for (modality in unique(data_siips_subregions_points_siips$modality)) {
      cur_data = data_siips_subregions_points_siips$score[data_siips_subregions_points_siips$stimulus_level == stimulus_level & data_siips_subregions_points_siips$condition == condition & data_siips_subregions_points_siips$modality == modality]
      range95 = quantile(cur_data, probs=c(.025, .975), na.rm = TRUE)
      data_siips_subregions_points_siips$min_range[data_siips_subregions_points_siips$stimulus_level == stimulus_level & data_siips_subregions_points_siips$condition == condition & data_siips_subregions_points_siips$modality == modality] = range95[1]
      data_siips_subregions_points_siips$max_range[data_siips_subregions_points_siips$stimulus_level == stimulus_level & data_siips_subregions_points_siips$condition == condition & data_siips_subregions_points_siips$modality == modality] = range95[2]
    }
  }
}
data_siips_subregions_points_siips_range95 = data_siips_subregions_points_siips[data_siips_subregions_points_siips$score >= data_siips_subregions_points_siips$min_range & data_siips_subregions_points_siips$score <= data_siips_subregions_points_siips$max_range,]
SIIPS_plot_95 = ggplot(data_for_plot_siips, aes(x=stimulus_level, y=mean, fill = condition)) + 
  geom_point(data=data_siips_subregions_points_siips_range95, aes(y = score, x=stimulus_level, color=condition), position=position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8), alpha=0.1, size = 0.8) +
  geom_bar(width=0.8,colour="black",position=position_dodge(0.8), stat="identity", alpha=0.1) + # Bar plot
  geom_errorbar(position=position_dodge(0.8), width=1/4, aes(ymin=mean-se_within, ymax=mean+se_within)) +
  geom_text(aes(y = max(mean) + 2 * sd, label=asterisks), size = 6, color = "black") +
  labs(x = "Stimulus level", y = "SIIPS score", color = "Condition", fill = "Condition") +
  facet_grid(. ~ modality, scales = "free", labeller = labeller(modality = modality.labs)) +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10)) +
  guides(color = "none", fill = guide_legend(title.position="top", title.hjust = 0.5))

## line plot, siips subregions
lineplot_siips_subregions = ggplot(data_siips_subregions_for_plot_no_siips, aes(x=stimulus_level, y=mean, colour=condition, group = condition)) + 
  geom_errorbar(aes(ymin=mean-se_within, ymax=mean+se_within), width=.1) +
  geom_line() +
  geom_point() +
  geom_text(aes(y = mean + se, label=asterisks), size = 6, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(x = "Stimulus level", y = "Mean score", colour = "Condition") +
  facet_grid(measure ~ modality, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10), legend.position="bottom")

lineplot_siips_subregions_1_4 = ggplot(data_siips_subregions_for_plot_no_siips1, aes(x=stimulus_level, y=mean, colour=condition, group = condition)) + 
  geom_errorbar(aes(ymin=mean-se_within, ymax=mean+se_within), width=.1) +
  geom_line() +
  geom_point() +
  geom_text(aes(y = mean + se, label=asterisks), size = 6, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(x = "Stimulus level", y = "Mean score", colour = "Condition") +
  facet_grid(measure ~ modality, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10), legend.position="bottom")

lineplot_siips_subregions_5_8 = ggplot(data_siips_subregions_for_plot_no_siips2, aes(x=stimulus_level, y=mean, colour=condition, group = condition)) + 
  geom_errorbar(aes(ymin=mean-se_within, ymax=mean+se_within), width=.1) +
  geom_line() +
  geom_point() +
  geom_text(aes(y = mean + se, label=asterisks), size = 6, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(x = "Stimulus level", y = "Mean score", colour = "Condition") +
  facet_grid(measure ~ modality, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10), legend.position="bottom")

lineplot_siips_subregions_combined = ggarrange(lineplot_siips_subregions_1_4, lineplot_siips_subregions_5_8,
                                               labels = c("", ""),
                                               ncol = 2, nrow = 1, 
                                               common.legend = TRUE, legend = "bottom")

# bar plot with dots - siips subregions
bar_dots_siips_subregions_plot = ggplot(data_siips_subregions_for_plot_no_siips, aes(x=stimulus_level, y=mean, fill = condition)) + 
  geom_bar(width=0.8,colour="black",position=position_dodge(0.8), stat="identity", alpha=0.1) + # Bar plot
  geom_errorbar(position=position_dodge(0.8), width=1/4, aes(ymin=mean-se_within, ymax=mean+se_within)) +
  geom_point(data=data_siips_subregions_points_no_siips, aes(y = score, x=stimulus_level, color=condition), position=position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8), alpha=0.1, size = 0.8) +
  geom_text(aes(y = mean + 2.5*sd, label=asterisks), size = 6, color = "black") +
  labs(x = "Stimulus level", y = "Pain rating", fill = "Condition") +
  facet_grid(measure ~ modality, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10)) +
  guides(color = "none")

# bar plot with dots - ONLY 90% RANGE OF VALUES
for(measure in unique(data_siips_subregions_points_no_siips$measure)) {
  for (stimulus_level in unique(data_siips_subregions_points_no_siips$stimulus_level)) {
    for (condition in unique(data_siips_subregions_points_no_siips$condition)) {
      for (modality in unique(data_siips_subregions_points_no_siips$modality)) {
        cur_data = data_siips_subregions_points_no_siips$score[data_siips_subregions_points_no_siips$measure == measure & data_siips_subregions_points_no_siips$stimulus_level == stimulus_level & data_siips_subregions_points_no_siips$condition == condition & data_siips_subregions_points_no_siips$modality == modality]
        ranges = quantile(cur_data, probs=c(.025, .05, .95, .975), na.rm = TRUE)
        data_siips_subregions_points_no_siips$min_range[data_siips_subregions_points_no_siips$measure == measure & data_siips_subregions_points_no_siips$stimulus_level == stimulus_level & data_siips_subregions_points_no_siips$condition == condition & data_siips_subregions_points_no_siips$modality == modality] = ranges[1]
        data_siips_subregions_points_no_siips$range5per[data_siips_subregions_points_no_siips$measure == measure & data_siips_subregions_points_no_siips$stimulus_level == stimulus_level & data_siips_subregions_points_no_siips$condition == condition & data_siips_subregions_points_no_siips$modality == modality] = ranges[2]
        data_siips_subregions_points_no_siips$range95per[data_siips_subregions_points_no_siips$measure == measure & data_siips_subregions_points_no_siips$stimulus_level == stimulus_level & data_siips_subregions_points_no_siips$condition == condition & data_siips_subregions_points_no_siips$modality == modality] = ranges[3]
        data_siips_subregions_points_no_siips$max_range[data_siips_subregions_points_no_siips$measure == measure & data_siips_subregions_points_no_siips$stimulus_level == stimulus_level & data_siips_subregions_points_no_siips$condition == condition & data_siips_subregions_points_no_siips$modality == modality] = ranges[4]
      }
    }
  }
}
data_siips_subregions_points_no_siips_range90 = data_siips_subregions_points_no_siips[data_siips_subregions_points_no_siips$score >= data_siips_subregions_points_no_siips$range5per & data_siips_subregions_points_no_siips$score <= data_siips_subregions_points_no_siips$range95per,]
data_siips_subregions_points_no_siips_range90_1 = data_siips_subregions_points_no_siips_range90[data_siips_subregions_points_no_siips_range90$measure %in% siips_subregions.labs[2:5],]
data_siips_subregions_points_no_siips_range90_2 = data_siips_subregions_points_no_siips_range90[data_siips_subregions_points_no_siips_range90$measure %in% siips_subregions.labs[6:9],]

bar_dots_siips_subregions_plot_range90 = ggplot(data_siips_subregions_for_plot_no_siips, aes(x=stimulus_level, y=mean, fill = condition)) + 
  geom_point(data=data_siips_subregions_points_no_siips_range90, aes(y = score, x=stimulus_level, color=condition), position=position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8), alpha=0.1, size = 0.8) +
  geom_bar(width=0.8,colour="black",position=position_dodge(0.8), stat="identity", alpha=0.1) + # Bar plot
  geom_errorbar(position=position_dodge(0.8), width=1/4, aes(ymin=mean-se_within, ymax=mean+se_within)) +
  geom_text(aes(y = mean + 1.5*sd, label=asterisks), size = 6, color = "black") +
  labs(x = "Stimulus level", y = "Pain rating", fill = "Condition") +
  facet_grid(measure ~ modality, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10)) +
  guides(color = "none")

bar_dots_siips_subregions_plot_range90_1 = ggplot(data_siips_subregions_for_plot_no_siips1, aes(x=stimulus_level, y=mean, fill = condition)) + 
  geom_point(data=data_siips_subregions_points_no_siips_range90_1, aes(y = score, x=stimulus_level, color=condition), position=position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8), alpha=0.1, size = 0.8) +
  geom_bar(width=0.8,colour="black",position=position_dodge(0.8), stat="identity", alpha=0.1) + # Bar plot
  geom_errorbar(position=position_dodge(0.8), width=1/4, aes(ymin=mean-se_within, ymax=mean+se_within)) +
  geom_text(aes(y = mean + 1.5*sd, label=asterisks), size = 6, color = "black") +
  labs(x = "Stimulus level", y = "Score", fill = "Condition") +
  facet_grid(measure ~ modality, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10)) +
  guides(color = "none")

bar_dots_siips_subregions_plot_range90_2 = ggplot(data_siips_subregions_for_plot_no_siips2, aes(x=stimulus_level, y=mean, fill = condition)) + 
  geom_point(data=data_siips_subregions_points_no_siips_range90_2, aes(y = score, x=stimulus_level, color=condition), position=position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8), alpha=0.1, size = 0.8) +
  geom_bar(width=0.8,colour="black",position=position_dodge(0.8), stat="identity", alpha=0.1) + # Bar plot
  geom_errorbar(position=position_dodge(0.8), width=1/4, aes(ymin=mean-se_within, ymax=mean+se_within)) +
  geom_text(aes(y = mean + 1.5*sd, label=asterisks), size = 6, color = "black") +
  labs(x = "Stimulus level", y = "Score", fill = "Condition") +
  facet_grid(measure ~ modality, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10)) +
  guides(color = "none")

bar_dots_siips_subregions_plot_range90_combined = ggarrange(bar_dots_siips_subregions_plot_range90_1, bar_dots_siips_subregions_plot_range90_2,
                                                            labels = c("", ""),
                                                            ncol = 2, nrow = 1, 
                                                            common.legend = TRUE, legend = "bottom")

## higher level regions
increase_regions = c("canlab2018_mean_Ctx_8C_L", "canlab2018_mean_Ctx_46_L", "canlab2018_mean_Ctx_p9_46v_R", "canlab2018_mean_Ctx_a47r_R", "canlab2018_mean_Ctx_a10p_R", "canlab2018_mean_Ctx_p10p_R", "canlab2018_mean_Ctx_11l_L", "canlab2018_mean_Ctx_11l_R", "canlab2018_mean_V_Striatum_L","canlab2018_mean_V_Striatum_R")
increase_regions.labs = c("dlPFC L 1", "dlPFC L 2", "dlPFC R", "Lateral OFC", "Anterior OFC 1", "Anterior OFC 2", "Mid-lateral OFC L", "Mid-lateral OFC R", "NAc L", "NAc R")
names(increase_regions.labs) = increase_regions

# condition data for increase regions
data_increase_regions_long = gather(select(data, c(participant_ID, family_ID, stimLvl_labels, prodicaine_labels, heat_labels, starts_with("canlab2018"))), measure, score, starts_with("canlab2018") & !ends_with("Bstem_PAG"), factor_key=FALSE)
for (region_ind in 1:length(increase_regions)) {
  data_increase_regions_long$measure[data_increase_regions_long$measure == increase_regions[region_ind]] = increase_regions.labs[region_ind]
  stats_condition_z$measure[stats_condition_z$measure == increase_regions[region_ind]] = increase_regions.labs[region_ind]
}
data_increase_regions_long$measure = factor(data_increase_regions_long$measure)
data_increase_regions_for_plot = aggregate(data_increase_regions_long$score,list(measure = data_increase_regions_long$measure, modality = data_increase_regions_long$heat_labels, condition = data_increase_regions_long$prodicaine_labels, stimulus_level = data_increase_regions_long$stimLvl_labels, participant_ID = data_increase_regions_long$participant_ID), function(x) c(mean = mean(x, na.rm=TRUE)))
data_increase_regions_for_plot = aggregate(data_increase_regions_for_plot$x,list(measure = data_increase_regions_for_plot$measure, modality = data_increase_regions_for_plot$modality, condition = data_increase_regions_for_plot$condition, stimulus_level = data_increase_regions_for_plot$stimulus_level), function(x) c(mean = mean(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE), n = sum(!is.na(x))))
data_increase_regions_for_plot = cbind(data_increase_regions_for_plot[,1:ncol(data_increase_regions_for_plot)-1], as.data.frame(data_increase_regions_for_plot$x))
data_increase_regions_for_plot$se = data_increase_regions_for_plot$sd/sqrt(data_increase_regions_for_plot$n)
data_increase_regions_points = data_increase_regions_long
data_increase_regions_points = dplyr::rename(data_increase_regions_points, c(condition = prodicaine_labels, modality = heat_labels, stimulus_level = stimLvl_labels))

# compute the adjusted within-subject variance for each region separately, otherwise the variance increases because of variance across regions
for (region_ind in 1:length(increase_regions)) {
  cur_values = summarySEwithin(data_increase_regions_long[data_increase_regions_long$measure == increase_regions.labs[region_ind],],  measurevar = "score", withinvars = c("stimLvl_labels", "prodicaine_labels", "heat_labels", "measure"), idvar = "participant_ID", na.rm=TRUE, conf.interval=.95)
  if (region_ind > 1) {
    data_for_plot_increase_regions_adjusted = rbind(data_for_plot_increase_regions_adjusted, cur_values)
  } else {
    data_for_plot_increase_regions_adjusted = cur_values
  }
}
data_for_plot_increase_regions_adjusted = dplyr::rename(data_for_plot_increase_regions_adjusted, c(condition = prodicaine_labels, modality = heat_labels, stimulus_level = stimLvl_labels, mean_within = score, se_within = se, sd_within = sd, ci_within = ci, N_within = N))
data_increase_regions_for_plot = merge(data_increase_regions_for_plot, data_for_plot_increase_regions_adjusted)

stats_condition_z_to_merge = stats_condition_z[stats_condition_z$measure %in% increase_regions.labs & stats_condition_z$effect == "placebo", c("measure", "modality", "asterisks")]
data_increase_regions_for_plot = merge(data_increase_regions_for_plot, stats_condition_z_to_merge)
data_increase_regions_for_plot$asterisks[data_increase_regions_for_plot$condition != "control" | data_increase_regions_for_plot$stimulus_level != "med"] = NA;

# bars with dots
bar_dots_increase_regions = ggplot(data_increase_regions_for_plot, aes(x=stimulus_level, y=mean, fill = condition)) + 
  geom_bar(width=0.8,colour="black",position=position_dodge(0.8), stat="identity", alpha=0.1) + # Bar plot
  geom_errorbar(position=position_dodge(0.8), width=1/4, aes(ymin=mean-se_within, ymax=mean+se_within)) +
  geom_point(data=data_increase_regions_points, aes(y = score, x=stimulus_level, color=condition), position=position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8), alpha=0.1, size = 0.8) +
  geom_text(aes(y = mean + 2.5*sd, label=asterisks), size = 6, color = "black") +
  labs(x = "Stimulus level", y = "Score", fill = "Condition") +
  facet_grid(measure ~ modality, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10)) +
  guides(color = "none")

# bar plot with dots - ONLY 90% RANGE OF VALUES
for(measure in unique(data_increase_regions_points$measure)) {
  for (stimulus_level in unique(data_increase_regions_points$stimulus_level)) {
    for (condition in unique(data_increase_regions_points$condition)) {
      for (modality in unique(data_increase_regions_points$modality)) {
        cur_data = data_increase_regions_points$score[data_increase_regions_points$measure == measure & data_increase_regions_points$stimulus_level == stimulus_level & data_increase_regions_points$condition == condition & data_increase_regions_points$modality == modality]
        ranges = quantile(cur_data, probs=c(.025, .05, .95, .975), na.rm = TRUE)
        data_increase_regions_points$min_range[data_increase_regions_points$measure == measure & data_increase_regions_points$stimulus_level == stimulus_level & data_increase_regions_points$condition == condition & data_increase_regions_points$modality == modality] = ranges[1]
        data_increase_regions_points$range5per[data_increase_regions_points$measure == measure & data_increase_regions_points$stimulus_level == stimulus_level & data_increase_regions_points$condition == condition & data_increase_regions_points$modality == modality] = ranges[2]
        data_increase_regions_points$range95per[data_increase_regions_points$measure == measure & data_increase_regions_points$stimulus_level == stimulus_level & data_increase_regions_points$condition == condition & data_increase_regions_points$modality == modality] = ranges[3]
        data_increase_regions_points$max_range[data_increase_regions_points$measure == measure & data_increase_regions_points$stimulus_level == stimulus_level & data_increase_regions_points$condition == condition & data_increase_regions_points$modality == modality] = ranges[4]
      }
    }
  }
}
data_increase_regions_points_range90 = data_increase_regions_points[data_increase_regions_points$score >= data_increase_regions_points$range5per & data_increase_regions_points$score <= data_increase_regions_points$range95per,]
data_increase_regions_points_range90_1 = data_increase_regions_points_range90[data_increase_regions_points_range90$measure %in% increase_regions.labs[1:ceil(length(increase_regions)/2)],]
data_increase_regions_points_range90_2 = data_increase_regions_points_range90[data_increase_regions_points_range90$measure %in% increase_regions.labs[ceil(length(increase_regions)/2)+1:length(increase_regions)],]

data_increase_regions_for_plot1 = data_increase_regions_for_plot[data_increase_regions_for_plot$measure %in% increase_regions.labs[1:ceil(length(increase_regions)/2)],]
data_increase_regions_for_plot2 = data_increase_regions_for_plot[data_increase_regions_for_plot$measure %in% increase_regions.labs[ceil(length(increase_regions)/2)+1:length(increase_regions)],]

bar_dots_increase_regions_range90 = ggplot(data_increase_regions_for_plot, aes(x=stimulus_level, y=mean, fill = condition)) + 
  geom_point(data=data_increase_regions_points_range90, aes(y = score, x=stimulus_level, color=condition), position=position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8), alpha=0.1, size = 0.8) +
  geom_bar(width=0.8,colour="black",position=position_dodge(0.8), stat="identity", alpha=0.1) + # Bar plot
  geom_errorbar(position=position_dodge(0.8), width=1/4, aes(ymin=mean-se_within, ymax=mean+se_within)) +
  geom_text(aes(y = mean + sd, label=asterisks), size = 6, color = "black") +
  labs(x = "Stimulus level", y = "Score", fill = "Condition") +
  facet_grid(measure ~ modality, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10)) +
  guides(color = "none")

bar_dots_increase_regions_range90_1 = ggplot(data_increase_regions_for_plot1, aes(x=stimulus_level, y=mean, fill = condition)) + 
  geom_point(data=data_increase_regions_points_range90_1, aes(y = score, x=stimulus_level, color=condition), position=position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8), alpha=0.1, size = 0.8) +
  geom_bar(width=0.8,colour="black",position=position_dodge(0.8), stat="identity", alpha=0.1) + # Bar plot
  geom_errorbar(position=position_dodge(0.8), width=1/4, aes(ymin=mean-se_within, ymax=mean+se_within)) +
  geom_text(aes(y = mean + 1.5*sd, label=asterisks), size = 6, color = "black") +
  labs(x = "Stimulus level", y = "Score", fill = "Condition") +
  facet_grid(measure ~ modality, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10)) +
  guides(color = "none")

bar_dots_increase_regions_range90_2 = ggplot(data_increase_regions_for_plot2, aes(x=stimulus_level, y=mean, fill = condition)) + 
  geom_point(data=data_increase_regions_points_range90_2, aes(y = score, x=stimulus_level, color=condition), position=position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8), alpha=0.1, size = 0.8) +
  geom_bar(width=0.8,colour="black",position=position_dodge(0.8), stat="identity", alpha=0.1) + # Bar plot
  geom_errorbar(position=position_dodge(0.8), width=1/4, aes(ymin=mean-se_within, ymax=mean+se_within)) +
  geom_text(aes(y = mean + 1.5*sd, label=asterisks), size = 6, color = "black") +
  labs(x = "Stimulus level", y = "Score", fill = "Condition") +
  facet_grid(measure ~ modality, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10)) +
  guides(color = "none")

bar_dots_increase_regions_range90_combined = ggarrange(bar_dots_increase_regions_range90_1, bar_dots_increase_regions_range90_2,
                                                       labels = c("", ""),
                                                       ncol = 2, nrow = 1, 
                                                       common.legend = TRUE, legend = "bottom")

# line plots, increased regions
lineplot_increase_regions = ggplot(data_increase_regions_for_plot, aes(x=stimulus_level, y=mean, colour=condition, group = condition)) + 
  geom_errorbar(aes(ymin=mean-se_within, ymax=mean+se_within), width=.1) +
  geom_line() +
  geom_point() +
  geom_text(aes(y = max(mean) + 0.02, label=asterisks), size = 6, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_y_continuous(breaks = seq(-0.2, 0.2, 0.1), limits = c(-0.2, 0.2)) +
  labs(x = "Stimulus level", y = "Mean score", colour = "Condition") +
  facet_grid(modality ~ measure, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 9), legend.position="bottom")

### create summarized df for plotting the correlation between behavioral and neural placebo-induced reductions
data_analgesia_long = gather(select(data_analgesia, c(participant_ID, family_ID, stimLvl_labels, heat_labels, Yint, nps, siips)), measure, score, c(nps, siips), factor_key=TRUE)

# NPS
# New facet label names for modality variable
modality.labs = c("Thermal", "Mechanical")
names(modality.labs) = c("thermal", "mechanical")
corr_plot_asterisks = select(stats_analgesia_z,c("measure", "modality", "asterisks"))
corr_plot_asterisks = dplyr::rename(corr_plot_asterisks, heat_labels = modality)
data_analgesia_long_all = gather(select(data_analgesia, c(participant_ID, family_ID, stimLvl_labels, heat_labels, Yint, nps:rank_siips)), measure, score, nps:rank_siips, factor_key=TRUE)
median_Yint_data = aggregate(data_analgesia_long_all$Yint, list(heat_labels = data_analgesia_long_all$heat_labels, measure = data_analgesia_long_all$measure), function(x) c(x = median(x, na.rm=TRUE)))
median_Yint_data = median_Yint_data[median_Yint_data$measure %in% corr_plot_asterisks$measure,]
median_Yint_data = dplyr::rename(median_Yint_data, median_Yint = x)
max_score_data = aggregate(data_analgesia_long_all$score, list(heat_labels = data_analgesia_long_all$heat_labels, measure = data_analgesia_long_all$measure), function(x) c(x = max(x, na.rm=TRUE)))
max_score_data = max_score_data[max_score_data$measure %in% corr_plot_asterisks$measure,]
max_score_data = dplyr::rename(max_score_data, max_score = x)
corr_plot_asterisks = merge(corr_plot_asterisks, median_Yint_data)     
corr_plot_asterisks = merge(corr_plot_asterisks, max_score_data)  

NPS_corr_plot = ggplot(data = data_analgesia_long[data_analgesia_long$measure == "nps",]) +
  geom_point(aes(x = Yint, y = score, colour = stimLvl_labels), size = 0.5, alpha = 0.1) +
  geom_smooth (aes(x = Yint, y = score, colour = stimLvl_labels), method = "lm", alpha=0.3, size=0) +
  stat_smooth (aes(x = Yint, y = score, colour = stimLvl_labels), geom="line", method = "lm", size=1, alpha=0.5) +
  geom_text(data = corr_plot_asterisks[corr_plot_asterisks$measure == "nps",], mapping = aes(x = median_Yint, y = max_score+1, label = asterisks)) +
  labs(x = "Pain ratings\nControl - Placebo", y = "NPS score\nControl - Placebo", colour = "Stimulus level") +
  scale_colour_manual(values = c("orange", "red2", "brown")) +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10)) +
  facet_grid(. ~ heat_labels, scales = "free_y", labeller = labeller(heat_labels = modality.labs)) +
  guides(color = guide_legend(title.position="top", title.hjust = 0.5))

# SIIPS
SIIPS_corr_plot = ggplot(data = data_analgesia_long[data_analgesia_long$measure == "siips",]) +
  geom_point(aes(x = Yint, y = score, colour = stimLvl_labels), size = 0.5, alpha = 0.1) +
  geom_smooth (aes(x = Yint, y = score, colour = stimLvl_labels), method = "lm", alpha=0.3, size=0) +
  stat_smooth (aes(x = Yint, y = score, colour = stimLvl_labels), geom="line", method = "lm", size=1, alpha=0.5) +
  geom_text(data = corr_plot_asterisks[corr_plot_asterisks$measure == "siips",], mapping = aes(x = median_Yint, y = max_score+1, label = asterisks)) +
  labs(x = "Pain ratings\nControl - Placebo", y = "SIIPS score\nControl - Placebo", colour = "Stimulus level") +
  scale_colour_manual(values = c("orange", "red2", "brown")) +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10)) +
  facet_grid(. ~ heat_labels, labeller = labeller(heat_labels = modality.labs), scales = "free") +
  guides(color = guide_legend(title.position="top", title.hjust = 0.5))

## correlation plots for all nociceptive regions
data_analgesia_long_nociceptive = gather(select(data_analgesia, c(participant_ID, family_ID, stimLvl_labels, heat_labels, Yint, starts_with("pain_pathways") | starts_with("canlab2018_mean_Bstem_PAG"))), measure, score, c(starts_with("pain_pathways") | starts_with("canlab2018_mean_Bstem_PAG")), factor_key=TRUE)
# New facet label names for nociceptive region variable
region_nociceptive.labs = c("aMCC", "dpIns L", "dpIns R", "PAG", "Med Thal", "VPL/M Thal L", "VPL/M Thal R")
names(region_nociceptive.labs) = c("pain_pathways_aMCC_MPFC", "pain_pathways_dpIns_L", "pain_pathways_dpIns_R", "canlab2018_mean_Bstem_PAG", "pain_pathways_Thal_MD", "pain_pathways_Thal_VPLM_L", "pain_pathways_Thal_VPLM_R")
data_analgesia_long_nociceptive$measure = factor(data_analgesia_long_nociceptive$measure, levels=c("pain_pathways_aMCC_MPFC", "pain_pathways_dpIns_L", "pain_pathways_dpIns_R", "canlab2018_mean_Bstem_PAG", "pain_pathways_Thal_MD", "pain_pathways_Thal_VPLM_L", "pain_pathways_Thal_VPLM_R"))
nociceptive_regions_corr_plot = ggplot(data = data_analgesia_long_nociceptive) +
  geom_point(size = 0.5, alpha = 0.1, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  geom_smooth (method = "lm", alpha=0.2, size=0, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  stat_smooth (geom="line", method = "lm", size=1, alpha=0.5, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  geom_text(data = corr_plot_asterisks[corr_plot_asterisks$measure %in% names(region_nociceptive.labs),], mapping = aes(x = median_Yint, y = max_score, label = asterisks)) +
  labs(x = "Pain ratings\nControl - Placebo", y = "Neural activity\nControl - Placebo", colour = "Stimulus level") +
  scale_colour_manual(values = c("orange", "red2", "brown")) +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10), legend.position="right", strip.text.y = element_text(angle = 0)) +
  facet_grid(measure ~ heat_labels, labeller = labeller(measure = region_nociceptive.labs), scales = "free")

nociceptive_regions_corr_plot_1 = ggplot(data = data_analgesia_long_nociceptive[data_analgesia_long_nociceptive$measure %in% levels(data_analgesia_long_nociceptive$measure)[1:ceil(length(nociceptive_regions)/2)],]) +
  geom_point(size = 0.5, alpha = 0.1, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  geom_smooth (method = "lm", alpha=0.2, size=0, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  stat_smooth (geom="line", method = "lm", size=1, alpha=0.5, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  geom_text(data = corr_plot_asterisks[corr_plot_asterisks$measure %in% levels(data_analgesia_long_nociceptive$measure)[1:ceil(length(nociceptive_regions)/2)],], mapping = aes(x = median_Yint, y = max_score, label = asterisks)) +
  labs(x = "Pain ratings\nControl - Placebo", y = "Neural activity\nControl - Placebo", colour = "Stimulus level") +
  scale_colour_manual(values = c("orange", "red2", "brown")) +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10), legend.position="right", strip.text.y = element_text(angle = 0)) +
  facet_grid(measure ~ heat_labels, labeller = labeller(measure = region_nociceptive.labs[1:ceil(length(nociceptive_regions)/2)]), scales = "free")

nociceptive_regions_corr_plot_2 = ggplot(data = data_analgesia_long_nociceptive[data_analgesia_long_nociceptive$measure %in% levels(data_analgesia_long_nociceptive$measure)[ceil(length(nociceptive_regions)/2)+1:length(nociceptive_regions)],]) +
  geom_point(size = 0.5, alpha = 0.1, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  geom_smooth (method = "lm", alpha=0.2, size=0, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  stat_smooth (geom="line", method = "lm", size=1, alpha=0.5, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  geom_text(data = corr_plot_asterisks[corr_plot_asterisks$measure %in% levels(data_analgesia_long_nociceptive$measure)[ceil(length(nociceptive_regions)/2)+1:length(nociceptive_regions)],], mapping = aes(x = median_Yint, y = max_score, label = asterisks)) +
  labs(x = "Pain ratings\nControl - Placebo", y = "Neural activity\nControl - Placebo", colour = "Stimulus level") +
  scale_colour_manual(values = c("orange", "red2", "brown")) +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10), legend.position="right", strip.text.y = element_text(angle = 0)) +
  facet_grid(measure ~ heat_labels, labeller = labeller(measure = region_nociceptive.labs[ceil(length(nociceptive_regions)/2)+1:length(nociceptive_regions)]), scales = "free")

# adding the following lines because the number of nociceptive regions is odd
nociceptive_regions_corr_plot_2_with_space = ggarrange(nociceptive_regions_corr_plot_2, NULL,
                                                       nrow = 2, heights = c(0.9,0.25),
                                                       common.legend = TRUE, legend = "none")
nociceptive_regions_corr_plot_combined = ggarrange(nociceptive_regions_corr_plot_1, nociceptive_regions_corr_plot_2_with_space,
                                                   labels = "", ncol = 2, nrow = 1,
                                                   common.legend = TRUE, legend = "bottom")

## correlation plots for siips sub-regions
data_analgesia_no_coss = select(data_analgesia, !ends_with("_coss"))
data_analgesia_long_siips_subregions = gather(select(data_analgesia_no_coss, c(participant_ID, family_ID, stimLvl_labels, heat_labels, Yint, siips, starts_with("siipspos"), starts_with("siipsneg"))), measure, score, c(siips, starts_with("siipspos"), starts_with("siipsneg")), factor_key=TRUE)
data_analgesia_long_siips_subregions = data_analgesia_long_siips_subregions[data_analgesia_long_siips_subregions$measure %in% siips_subregions,]
data_analgesia_long_siips_subregions$measure = factor(data_analgesia_long_siips_subregions$measure, levels=siips_subregions)
siips_subregions_corr_plot = ggplot(data = data_analgesia_long_siips_subregions) +
  geom_point(size = 0.5, alpha = 0.1, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  geom_smooth (method = "lm", alpha=0.2, size=0, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  stat_smooth (geom="line", method = "lm", size=1, alpha=0.4, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  geom_text(data = corr_plot_asterisks[corr_plot_asterisks$measure %in% names(siips_subregions.labs),], mapping = aes(x = median_Yint, y = max_score, label = asterisks)) +
  labs(x = "Pain ratings\nControl - Placebo", y = "Neural activity\nControl - Placebo", colour = "Stimulus level") +
  scale_colour_manual(values = c("orange", "red2", "brown")) +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10), legend.position="right", strip.text.y = element_text(angle = 0)) +
  facet_grid(measure ~ heat_labels, scales = "free" , labeller = labeller(measure = siips_subregions.labs))

siips_subregions_corr_plot_1 = ggplot(data = data_analgesia_long_siips_subregions[data_analgesia_long_siips_subregions$measure %in% levels(data_analgesia_long_siips_subregions$measure)[2:5],]) +
  geom_point(size = 0.5, alpha = 0.1, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  geom_smooth (method = "lm", alpha=0.2, size=0, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  stat_smooth (geom="line", method = "lm", size=1, alpha=0.4, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  geom_text(data = corr_plot_asterisks[corr_plot_asterisks$measure %in% levels(data_analgesia_long_siips_subregions$measure)[2:5],], mapping = aes(x = median_Yint, y = max_score, label = asterisks)) +
  labs(x = "Pain ratings\nControl - Placebo", y = "Neural activity\nControl - Placebo", colour = "Stimulus level") +
  scale_colour_manual(values = c("orange", "red2", "brown")) +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10), legend.position="right", strip.text.y = element_text(angle = 0)) +
  facet_grid(measure ~ heat_labels, scales = "free" , labeller = labeller(measure = siips_subregions.labs[2:5]))

siips_subregions_corr_plot_2 = ggplot(data = data_analgesia_long_siips_subregions[data_analgesia_long_siips_subregions$measure %in% levels(data_analgesia_long_siips_subregions$measure)[6:9],]) +
  geom_point(size = 0.5, alpha = 0.1, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  geom_smooth (method = "lm", alpha=0.2, size=0, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  stat_smooth (geom="line", method = "lm", size=1, alpha=0.4, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  geom_text(data = corr_plot_asterisks[corr_plot_asterisks$measure %in% levels(data_analgesia_long_siips_subregions$measure)[6:9],], mapping = aes(x = median_Yint, y = max_score, label = asterisks)) +
  labs(x = "Pain ratings\nControl - Placebo", y = "Neural activity\nControl - Placebo", colour = "Stimulus level") +
  scale_colour_manual(values = c("orange", "red2", "brown")) +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10), legend.position="right", strip.text.y = element_text(angle = 0)) +
  facet_grid(measure ~ heat_labels, scales = "free" , labeller = labeller(measure = siips_subregions.labs[6:9]))

siips_subregions_corr_plot_combined = ggarrange(siips_subregions_corr_plot_1, siips_subregions_corr_plot_2,
                                                labels = "", ncol = 2, nrow = 1,
                                                common.legend = TRUE, legend = "bottom")

## correlation plots for highr level ROIs
data_analgesia_long_increase_regions = gather(select(data_analgesia, c(participant_ID, family_ID, stimLvl_labels, heat_labels, Yint, starts_with("canlab2018") & !ends_with("Bstem_PAG"))), measure, score, c(starts_with("canlab2018")), factor_key=TRUE)
data_analgesia_long_increase_regions$measure = factor(data_analgesia_long_increase_regions$measure, levels=increase_regions)
increase_regions_corr_plot = ggplot(data = data_analgesia_long_increase_regions) +
  geom_point(size = 0.5, alpha = 0.1, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  geom_smooth (method = "lm", alpha=0.2, size=0, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  stat_smooth (geom="line", method = "lm", size=1, alpha=0.4, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  geom_text(data = corr_plot_asterisks[corr_plot_asterisks$measure %in% names(increase_regions.labs),], mapping = aes(x = median_Yint, y = max_score, label = asterisks)) +
  labs(x = "Pain ratings\nControl - Placebo", y = "Neural activity\nControl - Placebo", colour = "Stimulus level") +
  scale_colour_manual(values = c("orange", "red2", "brown")) +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10), legend.position="right", strip.text.y = element_text(angle = 0)) +
  facet_grid(measure ~ heat_labels, scales = "free" , labeller = labeller(measure = increase_regions.labs))

increase_regions1 = increase_regions[1:ceil(length(increase_regions)/2)]
increase_regions.labs1 = increase_regions.labs[1:ceil(length(increase_regions.labs)/2)]
increase_regions_corr_plot_1 = ggplot(data = data_analgesia_long_increase_regions[data_analgesia_long_increase_regions$measure %in% increase_regions1,]) +
  geom_point(size = 0.5, alpha = 0.1, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  geom_smooth (method = "lm", alpha=0.2, size=0, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  stat_smooth (geom="line", method = "lm", size=1, alpha=0.4, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  geom_text(data = corr_plot_asterisks[corr_plot_asterisks$measure %in% increase_regions1,], mapping = aes(x = median_Yint, y = max_score, label = asterisks)) +
  labs(x = "Pain ratings\nControl - Placebo", y = "Neural activity\nControl - Placebo", colour = "Stimulus level") +
  scale_colour_manual(values = c("orange", "red2", "brown")) +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10), legend.position="right", strip.text.y = element_text(angle = 0)) +
  facet_grid(measure ~ heat_labels, scales = "free_y" , labeller = labeller(measure = increase_regions.labs1))

increase_regions2 = increase_regions[ceil(length(increase_regions)/2)+1:length(increase_regions)]
increase_regions.labs2 = increase_regions.labs[ceil(length(increase_regions.labs)/2)+1:length(increase_regions.labs)]
increase_regions_corr_plot_2 = ggplot(data = data_analgesia_long_increase_regions[data_analgesia_long_increase_regions$measure %in% increase_regions2,]) +
  geom_point(size = 0.5, alpha = 0.1, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  geom_smooth (method = "lm", alpha=0.2, size=0, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  stat_smooth (geom="line", method = "lm", size=1, alpha=0.4, aes(x = Yint, y = score, colour = stimLvl_labels)) +
  geom_text(data = corr_plot_asterisks[corr_plot_asterisks$measure %in% increase_regions2,], mapping = aes(x = median_Yint, y = max_score, label = asterisks)) +
  labs(x = "Pain ratings\nControl - Placebo", y = "Neural activity\nControl - Placebo", colour = "Stimulus level") +
  scale_colour_manual(values = c("orange", "red2", "brown")) +
  theme_bw() +
  theme(strip.background = element_rect(color="white", fill="white"), text = element_text(size = 10), legend.position="right", strip.text.y = element_text(angle = 0)) +
  facet_grid(measure ~ heat_labels, scales = "free_y" , labeller = labeller(measure = increase_regions.labs2))

increase_regions_corr_plot_combined = ggarrange(increase_regions_corr_plot_1, increase_regions_corr_plot_2,
                                                labels = "", ncol = 2, nrow = 1,
                                                common.legend = TRUE, legend = "bottom")


#######################
### prepare figures ###
#######################

## figure 2 - behavioral results
figure2 = lineplot_behav

# save figure
size_factor = 1
ggsave('fig2.pdf',
       plot = figure2,
       path = fig_output_path,
       width = 120*size_factor,
       height = 100*size_factor,
       units = "mm",
       dpi = 300)

## figure 3 - NPS + nociceptive regions
# combine the nociceptive plots together
NPS_brain = readPNG("figures/NPS_surface_contrast30_slices_contrast10.png")
gg_NPS_brain = ggplot() + 
  background_image(NPS_brain) +
  # This ensures that the image leaves some space at the edges
  theme(plot.margin = margin(t=0.2, l=0.2, r=0.2, b=0.2, unit = "cm"), panel.background = element_rect(fill = "white"))

NPS_brain_and_plots = ggarrange(gg_NPS_brain, lineplot_nps_ver,
                                labels = c("A", "B"),
                                ncol = 2, nrow = 1,
                                widths = c(2,1))

figure3 = ggarrange(NPS_brain_and_plots, lineplot_nociceptive_regions,
                    labels = c("" , "C"),
                    ncol = 1, nrow = 2,
                    #heights = c(1.2, 2),
                    common.legend = FALSE, legend = "bottom")

# save figure
size_factor = 1
ggsave('fig3.pdf',
       plot = figure3,
       path = fig_output_path,
       width = 180*size_factor,
       height = 210*size_factor,
       units = "mm",
       dpi = 300)

## figure 4 - siips
SIIPS_brain = readPNG("figures/SIIPS_multi_surface_contrast12_5_slices_contrast5.png")
gg_SIIPS_brain = ggplot() + 
  background_image(SIIPS_brain) +
  # This ensures that the image leaves some space at the edges
  theme(plot.margin = margin(t=0.1, l=0.1, r=0.1, b=0.1, unit = "cm"), panel.background = element_rect(fill = "white"))

SIIPS_brain_and_plots = ggarrange(gg_SIIPS_brain, lineplot_siips_ver,
                                  labels = c("A", "B"),
                                  ncol = 2, nrow = 1,
                                  widths = c(2.5,1))

figure4 = SIIPS_brain_and_plots

# save figure
size_factor = 1
ggsave('fig4.pdf',
       plot = figure4,
       path = fig_output_path,
       width = 210*size_factor,
       height = 100*size_factor,
       units = "mm",
       dpi = 300)

## Figure 5 - high-level regions
figure5 = lineplot_increase_regions

# save figure
size_factor = 1
ggsave('fig5.pdf',
       plot = figure5,
       path = fig_output_path,
       width = 240*size_factor,
       height = 100*size_factor,
       units = "mm",
       dpi = 300)

## Figure 6 - NPS and SIIPS corr plots
nps_siips_corr = ggarrange(NPS_corr_plot, SIIPS_corr_plot,
                           labels = c("A", "B"),
                           ncol = 2, nrow = 1,
                           common.legend = TRUE, legend = "bottom")

figure6 = nps_siips_corr

# save figure
size_factor = 1
ggsave('fig6.pdf',
       plot = figure6,
       path = fig_output_path,
       width = 240*size_factor,
       height = 80*size_factor,
       units = "mm",
       dpi = 300)



## Supp Figures
# Supp Figure 1
supp_fig1 = bar_dots_behav_effect_plot
size_factor = 1
ggsave('supp_fig1.pdf',
       plot = supp_fig1,
       path = fig_output_path,
       width = 120*size_factor,
       height = 80*size_factor,
       units = "mm",
       dpi = 300)

warning("supp figure 2- NPS-\nthe dots included are only from the 95% mid-range of each category- 2.5% top and 2.5% bottom values are excluded from visualization.\nThe bars are based on the mean of *all* values, 100% of the range")
supp_fig2 = NPS_plot_95
size_factor = 1
ggsave('supp_fig2.pdf',
       plot = supp_fig2,
       path = fig_output_path,
       width = 130*size_factor,
       height = 100*size_factor,
       units = "mm",
       dpi = 300)

warning("supp figure 3- nociceptive regions-\nthe dots included are only from the 90% mid-range of each category- 5% top and 5% bottom values are excluded from visualization.\nThe bars are based on the mean of *all* values, 100% of the range")
supp_fig3 = bar_dots_nociceptive_regions_plot_range90_combined
size_factor = 1
ggsave('supp_fig3.pdf',
       plot = supp_fig3,
       path = fig_output_path,
       width = 160*size_factor,
       height = 210*size_factor,
       units = "mm",
       dpi = 300)

warning("supp figure 4- SIIPS-\nthe dots included are only from the 95% mid-range of each category- 2.5% top and 2.5% bottom values are excluded from visualization.\nThe bars are based on the mean of *all* values, 100% of the range")
supp_fig4 = SIIPS_plot_95
size_factor = 1
ggsave('supp_fig4.pdf',
       plot = supp_fig4,
       path = fig_output_path,
       width = 130*size_factor,
       height = 100*size_factor,
       units = "mm",
       dpi = 300)

supp_fig5 = lineplot_siips_subregions_combined
size_factor = 1
ggsave('supp_fig5.pdf',
       plot = supp_fig5,
       path = fig_output_path,
       width = 180*size_factor,
       height = 200*size_factor,
       units = "mm",
       dpi = 300)

warning("supp figure 6- SIIPS subregions-\nthe dots included are only from the 90% mid-range of each category- 5% top and 5% bottom values are excluded from visualization.\nThe bars are based on the mean of *all* values, 100% of the range")
supp_fig6 = bar_dots_siips_subregions_plot_range90_combined
size_factor = 1
ggsave('supp_fig6.pdf',
       plot = supp_fig6,
       path = fig_output_path,
       width = 160*size_factor,
       height = 210*size_factor,
       units = "mm",
       dpi = 300)

warning("supp figure 7- higher level regions-\nthe dots included are only from the 90% mid-range of each category- 5% top and 5% bottom values are excluded from visualization.\nThe bars are based on the mean of *all* values, 100% of the range")
supp_fig7 = bar_dots_increase_regions_range90_combined
size_factor = 1
ggsave('supp_fig7.pdf',
       plot = supp_fig7,
       path = fig_output_path,
       width = 160*size_factor,
       height = 210*size_factor,
       units = "mm",
       dpi = 300)



#############
# GET MEANS #
#############
# summarize data to get means of pain ratings in each combination of modality and stimulus level
Yint_means = aggregate(data$Yint,list(modality = data$heat_labels, condition = data$prodicaine_labels, stimulus_level = data$stimLvl_labels, participant_ID = data$participant_ID), function(x) c(mean = mean(x, na.rm=TRUE)))
Yint_means = aggregate(Yint_means$x,list(modality = Yint_means$modality, condition = Yint_means$condition, stimulus_level = Yint_means$stimulus_level), function(x) c(mean = mean(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE), n = sum(!is.na(x))))
Yint_means = cbind(Yint_means[,1:ncol(Yint_means)-1], as.data.frame(Yint_means$x))
Yint_means$se = Yint_means$sd/sqrt(Yint_means$n)
Yint_means_stim_level = aggregate(Yint_means$mean, list(modality = Yint_means$modality, stimulus_level = Yint_means$stimulus_level), function(x) c(mean = mean(x, na.rm=TRUE)))
Yint_means_placebo_condition = aggregate(Yint_means$mean, list(modality = Yint_means$modality, condition = Yint_means$condition), function(x) c(mean = mean(x, na.rm=TRUE)))

# behavioral analgesia by sex
Yint_means_sex = aggregate(data_analgesia$Yint,list(modality = data_analgesia$heat_labels, sex = data_analgesia$sex, participant_ID = data_analgesia$participant_ID), function(x) c(mean = mean(x, na.rm=TRUE)))
Yint_means_sex = aggregate(Yint_means_sex$x,list(modality = Yint_means_sex$modality, sex = Yint_means_sex$sex), function(x) c(mean = mean(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE), n = sum(!is.na(x))))

# Test correlation between intensity and unpleasantness ratings
# first, without accounting for participants
cor.test(data$Yint, data$Yunp)
# compute the correlation for each participant
participant_IDs = unique(data$participant_ID)
cor_yint_yunp = as.data.frame(participant_IDs)
for (ind in 1:length(participant_IDs)) {
  cur_participant_ID = participant_IDs[ind]
  cur_data = data[data$participant_ID == cur_participant_ID, ]
  cur_data = cur_data[complete.cases(cur_data[ , c('Yint', 'Yunp')]), ]
  if (nrow(cur_data) > 1) {
    cor_yint_yunp$cor[ind] = cor(cur_data$Yint, cur_data$Yunp)
  } else {
    cor_yint_yunp$cor[ind] = NA
  }
}
mean_cor_yint_yunp = mean(cor_yint_yunp$cor, na.rm = T)
median_cor_yint_yunp = median(cor_yint_yunp$cor, na.rm = T)
sd_cor_yint_yunp = sd(cor_yint_yunp$cor, na.rm = T)
num_subjs_with_corr = sum(!is.na(cor_yint_yunp$cor))

## get means for unpleasantness ratings
Yunp_means = aggregate(data$Yunp,list(modality = data$heat_labels, condition = data$prodicaine_labels, stimulus_level = data$stimLvl_labels, participant_ID = data$participant_ID), function(x) c(mean = mean(x, na.rm=TRUE)))
Yunp_means = aggregate(Yunp_means$x,list(modality = Yunp_means$modality, condition = Yunp_means$condition, stimulus_level = Yunp_means$stimulus_level), function(x) c(mean = mean(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE), n = sum(!is.na(x))))
Yunp_means = cbind(Yunp_means[,1:ncol(Yunp_means)-1], as.data.frame(Yunp_means$x))
Yunp_means$se = Yunp_means$sd/sqrt(Yunp_means$n)
Yunp_means_stim_level = aggregate(Yunp_means$mean, list(modality = Yunp_means$modality, stimulus_level = Yunp_means$stimulus_level), function(x) c(mean = mean(x, na.rm=TRUE)))
Yunp_means_placebo_condition = aggregate(Yunp_means$mean, list(modality = Yunp_means$modality, condition = Yunp_means$condition), function(x) c(mean = mean(x, na.rm=TRUE)))