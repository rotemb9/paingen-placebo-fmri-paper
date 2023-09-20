# This script performs the Bayes Factor analysis
# Written by Rotem Botvinik-Nezer

library(BayesFactor)
library(dplyr)

# clear workspace
rm(list=ls())

# read the data
use_z_scored_data = FALSE
if (use_z_scored_data) {
  print("Using z scored data")
  data = read.csv("pain_evoked_brain_measures_z.csv")
} else {
  print("Using NON z scored data")
  data = read.csv("pain_evoked_brain_measures.csv")  
}
data$participant_ID = as.factor(data$participant_ID)
data$family_ID = as.factor(data$family_ID)
data$stimLvl = as.factor(data$stimLvl)
data$prodicaine = as.factor(data$prodicaine)
data$heat = as.factor(data$heat)

# run models for NPS
NPS_BF_test_heat = anovaBF(nps ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
NPS_BF_test_press = anovaBF(nps ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
NPS_BF_test_wide_heat = anovaBF(nps ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
NPS_BF_test_wide_press = anovaBF(nps ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
NPS_BF_test_ultrawide_heat = anovaBF(nps ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")
NPS_BF_test_ultrawide_press = anovaBF(nps ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")

# compare nps ~ prodicaine + stimLvl + participant_ID + family_ID to nps ~ stimLvl + participant_ID + family_ID (the null model for no effect of prodicaine on average, but random across participants)
NPS_BF_wide_no_int_heat = NPS_BF_test_wide_heat[3]/NPS_BF_test_wide_heat[1]
NPS_BF_wide_no_int_press = NPS_BF_test_wide_press[3]/NPS_BF_test_wide_press[1]
NPS_BF_ultrawide_no_int_heat = NPS_BF_test_ultrawide_heat[3]/NPS_BF_test_ultrawide_heat[1]
NPS_BF_ultrawide_no_int_press = NPS_BF_test_ultrawide_press[3]/NPS_BF_test_ultrawide_press[1]
NPS_BF_no_int_heat = NPS_BF_test_heat[3]/NPS_BF_test_heat[1]
NPS_BF_no_int_press = NPS_BF_test_press[3]/NPS_BF_test_press[1]
# compare nps ~ prodicaine * stimLvl + participant_ID + family_ID to nps ~ stimLvl + participant_ID + family_ID (the null model for no effect of prodicaine on average, but random across participants, with interaction in the alternative model)
NPS_BF_wide_int_heat = NPS_BF_test_wide_heat[4]/NPS_BF_test_wide_heat[1]
NPS_BF_wide_int_press = NPS_BF_test_wide_press[4]/NPS_BF_test_wide_press[1]
NPS_BF_ultrawide_int_heat = NPS_BF_test_ultrawide_heat[4]/NPS_BF_test_ultrawide_heat[1]
NPS_BF_ultrawide_int_press = NPS_BF_test_ultrawide_press[4]/NPS_BF_test_ultrawide_press[1]
NPS_BF_int_heat = NPS_BF_test_heat[4]/NPS_BF_test_heat[1]
NPS_BF_int_press = NPS_BF_test_press[4]/NPS_BF_test_press[1]

## SIIPS
siips_BF_test_heat = anovaBF(siips ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
siips_BF_test_press = anovaBF(siips ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
siips_BF_test_wide_heat = anovaBF(siips ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
siips_BF_test_wide_press = anovaBF(siips ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
siips_BF_test_ultrawide_heat = anovaBF(siips ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")
siips_BF_test_ultrawide_press = anovaBF(siips ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")

siips_BF_wide_no_int_heat = siips_BF_test_wide_heat[3]/siips_BF_test_wide_heat[1]
siips_BF_wide_no_int_press = siips_BF_test_wide_press[3]/siips_BF_test_wide_press[1]
siips_BF_wide_int_heat = siips_BF_test_wide_heat[4]/siips_BF_test_wide_heat[1]
siips_BF_wide_int_press = siips_BF_test_wide_press[4]/siips_BF_test_wide_press[1]
siips_BF_no_int_heat = siips_BF_test_heat[3]/siips_BF_test_heat[1]
siips_BF_no_int_press = siips_BF_test_press[3]/siips_BF_test_press[1]
siips_BF_ultrawide_no_int_heat = siips_BF_test_ultrawide_heat[3]/siips_BF_test_ultrawide_heat[1]
siips_BF_ultrawide_no_int_press = siips_BF_test_ultrawide_press[3]/siips_BF_test_ultrawide_press[1]

### a priori ROIs
# nociceptive (expected decrease)
# Dorsal posterior insula - left
dpIns_L_BF_test_heat = anovaBF(pain_pathways_dpIns_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
dpIns_L_BF_test_press = anovaBF(pain_pathways_dpIns_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
dpIns_L_BF_test_wide_heat = anovaBF(pain_pathways_dpIns_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
dpIns_L_BF_test_wide_press = anovaBF(pain_pathways_dpIns_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
dpIns_L_BF_test_ultrawide_heat = anovaBF(pain_pathways_dpIns_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")
dpIns_L_BF_test_ultrawide_press = anovaBF(pain_pathways_dpIns_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")

# Dorsal posterior insula - right
dpIns_R_BF_test_heat = anovaBF(pain_pathways_dpIns_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
dpIns_R_BF_test_press = anovaBF(pain_pathways_dpIns_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
dpIns_R_BF_test_wide_heat = anovaBF(pain_pathways_dpIns_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
dpIns_R_BF_test_wide_press = anovaBF(pain_pathways_dpIns_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
dpIns_R_BF_test_ultrawide_heat = anovaBF(pain_pathways_dpIns_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")
dpIns_R_BF_test_ultrawide_press = anovaBF(pain_pathways_dpIns_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")

# Anterior mid cingulate cortex
aMCC_MPFC_BF_test_heat = anovaBF(pain_pathways_aMCC_MPFC ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
aMCC_MPFC_BF_test_press = anovaBF(pain_pathways_aMCC_MPFC ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
aMCC_MPFC_BF_test_wide_heat = anovaBF(pain_pathways_aMCC_MPFC ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
aMCC_MPFC_BF_test_wide_press = anovaBF(pain_pathways_aMCC_MPFC ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
aMCC_MPFC_BF_test_ultrawide_heat = anovaBF(pain_pathways_aMCC_MPFC ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")
aMCC_MPFC_BF_test_ultrawide_press = anovaBF(pain_pathways_aMCC_MPFC ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")

# VPLM thalamus - left
Thal_VPLM_L_BF_test_heat = anovaBF(pain_pathways_Thal_VPLM_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
Thal_VPLM_L_BF_test_press = anovaBF(pain_pathways_Thal_VPLM_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
Thal_VPLM_L_BF_test_wide_heat = anovaBF(pain_pathways_Thal_VPLM_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
Thal_VPLM_L_BF_test_wide_press = anovaBF(pain_pathways_Thal_VPLM_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
Thal_VPLM_L_BF_test_ultrawide_heat = anovaBF(pain_pathways_Thal_VPLM_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")
Thal_VPLM_L_BF_test_ultrawide_press = anovaBF(pain_pathways_Thal_VPLM_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")

# VPLM thalamus - right
Thal_VPLM_R_BF_test_heat = anovaBF(pain_pathways_Thal_VPLM_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
Thal_VPLM_R_BF_test_press = anovaBF(pain_pathways_Thal_VPLM_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
Thal_VPLM_R_BF_test_wide_heat = anovaBF(pain_pathways_Thal_VPLM_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
Thal_VPLM_R_BF_test_wide_press = anovaBF(pain_pathways_Thal_VPLM_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
Thal_VPLM_R_BF_test_ultrawide_heat = anovaBF(pain_pathways_Thal_VPLM_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")
Thal_VPLM_R_BF_test_ultrawide_press = anovaBF(pain_pathways_Thal_VPLM_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")

# Medial thalamus
Thal_MD_BF_test_heat = anovaBF(pain_pathways_Thal_MD ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
Thal_MD_BF_test_press = anovaBF(pain_pathways_Thal_MD ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
Thal_MD_BF_test_wide_heat = anovaBF(pain_pathways_Thal_MD ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
Thal_MD_BF_test_wide_press = anovaBF(pain_pathways_Thal_MD ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
Thal_MD_BF_test_ultrawide_heat = anovaBF(pain_pathways_Thal_MD ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")
Thal_MD_BF_test_ultrawide_press = anovaBF(pain_pathways_Thal_MD ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")

# Dorsal posterior insula - left
PAG_BF_test_heat = anovaBF(canlab2018_mean_Bstem_PAG ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
PAG_BF_test_press = anovaBF(canlab2018_mean_Bstem_PAG ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
PAG_BF_test_wide_heat = anovaBF(canlab2018_mean_Bstem_PAG ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
PAG_BF_test_wide_press = anovaBF(canlab2018_mean_Bstem_PAG ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
PAG_BF_test_ultrawide_heat = anovaBF(canlab2018_mean_Bstem_PAG ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")
PAG_BF_test_ultrawide_press = anovaBF(canlab2018_mean_Bstem_PAG ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")

# BFs for a-priori nociceptive ROIs, wide, no interaction
dpIns_L_BF_wide_no_int_heat = dpIns_L_BF_test_wide_heat[3]/dpIns_L_BF_test_wide_heat[1]
dpIns_L_BF_wide_no_int_press = dpIns_L_BF_test_wide_press[3]/dpIns_L_BF_test_wide_press[1]
dpIns_R_BF_wide_no_int_heat = dpIns_R_BF_test_wide_heat[3]/dpIns_R_BF_test_wide_heat[1]
dpIns_R_BF_wide_no_int_press = dpIns_R_BF_test_wide_press[3]/dpIns_R_BF_test_wide_press[1]
aMCC_MPFC_BF_wide_no_int_heat = aMCC_MPFC_BF_test_wide_heat[3]/aMCC_MPFC_BF_test_wide_heat[1]
aMCC_MPFC_BF_wide_no_int_press = aMCC_MPFC_BF_test_wide_press[3]/aMCC_MPFC_BF_test_wide_press[1]
Thal_VPLM_L_BF_wide_no_int_heat = Thal_VPLM_L_BF_test_wide_heat[3]/Thal_VPLM_L_BF_test_wide_heat[1]
Thal_VPLM_L_BF_wide_no_int_press = Thal_VPLM_L_BF_test_wide_press[3]/Thal_VPLM_L_BF_test_wide_press[1]
Thal_VPLM_R_BF_wide_no_int_heat = Thal_VPLM_R_BF_test_wide_heat[3]/Thal_VPLM_R_BF_test_wide_heat[1]
Thal_VPLM_R_BF_wide_no_int_press = Thal_VPLM_R_BF_test_wide_press[3]/Thal_VPLM_R_BF_test_wide_press[1]
Thal_MD_BF_wide_no_int_heat = Thal_MD_BF_test_wide_heat[3]/Thal_MD_BF_test_wide_heat[1]
Thal_MD_BF_wide_no_int_press = Thal_MD_BF_test_wide_press[3]/Thal_MD_BF_test_wide_press[1]
PAG_BF_wide_no_int_heat = PAG_BF_test_wide_heat[3]/PAG_BF_test_wide_heat[1]
PAG_BF_wide_no_int_press = PAG_BF_test_wide_press[3]/PAG_BF_test_wide_press[1]

# put all in a list
list_nociception_bf_wide_no_int_heat = lapply(list(NPS = NPS_BF_wide_no_int_heat,
                                                    dpIns_L = dpIns_L_BF_wide_no_int_heat,
                                                    dpIns_R = dpIns_R_BF_wide_no_int_heat,
                                                    aMCC_MPFC = aMCC_MPFC_BF_wide_no_int_heat,
                                                    Thal_VPLM_L = Thal_VPLM_L_BF_wide_no_int_heat,
                                                    Thal_VPLM_R = Thal_VPLM_R_BF_wide_no_int_heat,
                                                    Thal_MD = Thal_MD_BF_wide_no_int_heat,
                                                    PAG = PAG_BF_wide_no_int_heat), extractBF)
list_nociception_bf_wide_no_int_press = lapply(list(NPS = NPS_BF_wide_no_int_press,
                                                    dpIns_L = dpIns_L_BF_wide_no_int_press,
                                                    dpIns_R = dpIns_R_BF_wide_no_int_press,
                                                    aMCC_MPFC = aMCC_MPFC_BF_wide_no_int_press,
                                                    Thal_VPLM_L = Thal_VPLM_L_BF_wide_no_int_press,
                                                    Thal_VPLM_R = Thal_VPLM_R_BF_wide_no_int_press,
                                                    Thal_MD = Thal_MD_BF_wide_no_int_press,
                                                    PAG = PAG_BF_wide_no_int_press), extractBF)

# robustness tests:
# with interaction
dpIns_L_BF_wide_int_heat = dpIns_L_BF_test_wide_heat[4]/dpIns_L_BF_test_wide_heat[1]
dpIns_L_BF_wide_int_press = dpIns_L_BF_test_wide_press[4]/dpIns_L_BF_test_wide_press[1]
dpIns_R_BF_wide_int_heat = dpIns_R_BF_test_wide_heat[4]/dpIns_R_BF_test_wide_heat[1]
dpIns_R_BF_wide_int_press = dpIns_R_BF_test_wide_press[4]/dpIns_R_BF_test_wide_press[1]
aMCC_MPFC_BF_wide_int_heat = aMCC_MPFC_BF_test_wide_heat[4]/aMCC_MPFC_BF_test_wide_heat[1]
aMCC_MPFC_BF_wide_int_press = aMCC_MPFC_BF_test_wide_press[4]/aMCC_MPFC_BF_test_wide_press[1]
Thal_VPLM_L_BF_wide_int_heat = Thal_VPLM_L_BF_test_wide_heat[4]/Thal_VPLM_L_BF_test_wide_heat[1]
Thal_VPLM_L_BF_wide_int_press = Thal_VPLM_L_BF_test_wide_press[4]/Thal_VPLM_L_BF_test_wide_press[1]
Thal_VPLM_R_BF_wide_int_heat = Thal_VPLM_R_BF_test_wide_heat[4]/Thal_VPLM_R_BF_test_wide_heat[1]
Thal_VPLM_R_BF_wide_int_press = Thal_VPLM_R_BF_test_wide_press[4]/Thal_VPLM_R_BF_test_wide_press[1]
Thal_MD_BF_wide_int_heat = Thal_MD_BF_test_wide_heat[4]/Thal_MD_BF_test_wide_heat[1]
Thal_MD_BF_wide_int_press = Thal_MD_BF_test_wide_press[4]/Thal_MD_BF_test_wide_press[1]
PAG_BF_wide_int_heat = PAG_BF_test_wide_heat[4]/PAG_BF_test_wide_heat[1]
PAG_BF_wide_int_press = PAG_BF_test_wide_press[4]/PAG_BF_test_wide_press[1]

# put all in a list
list_nociception_bf_wide_int_heat = lapply(list(NPS = NPS_BF_wide_int_heat,
                                                    dpIns_L = dpIns_L_BF_wide_int_heat,    
                                                    dpIns_R = dpIns_R_BF_wide_int_heat,
                                                    aMCC_MPFC = aMCC_MPFC_BF_wide_int_heat,
                                                    Thal_VPLM_L = Thal_VPLM_L_BF_wide_int_heat,
                                                    Thal_VPLM_R = Thal_VPLM_R_BF_wide_int_heat,
                                                    Thal_MD = Thal_MD_BF_wide_int_heat,
                                                    PAG = PAG_BF_wide_int_heat), extractBF)
list_nociception_bf_wide_int_press = lapply(list(NPS = NPS_BF_wide_int_press,
                                                     dpIns_L = dpIns_L_BF_wide_int_press,    
                                                     dpIns_R = dpIns_R_BF_wide_int_press,
                                                     aMCC_MPFC = aMCC_MPFC_BF_wide_int_press,
                                                     Thal_VPLM_L = Thal_VPLM_L_BF_wide_int_press,
                                                     Thal_VPLM_R = Thal_VPLM_R_BF_wide_int_press,
                                                     Thal_MD = Thal_MD_BF_wide_int_press,
                                                     PAG = PAG_BF_wide_int_press), extractBF)

# default width:
dpIns_L_BF_no_int_heat = dpIns_L_BF_test_heat[3]/dpIns_L_BF_test_heat[1]
dpIns_L_BF_no_int_press = dpIns_L_BF_test_press[3]/dpIns_L_BF_test_press[1]
dpIns_R_BF_no_int_heat = dpIns_R_BF_test_heat[3]/dpIns_R_BF_test_heat[1]
dpIns_R_BF_no_int_press = dpIns_R_BF_test_press[3]/dpIns_R_BF_test_press[1]
aMCC_MPFC_BF_no_int_heat = aMCC_MPFC_BF_test_heat[3]/aMCC_MPFC_BF_test_heat[1]
aMCC_MPFC_BF_no_int_press = aMCC_MPFC_BF_test_press[3]/aMCC_MPFC_BF_test_press[1]
Thal_VPLM_L_BF_no_int_heat = Thal_VPLM_L_BF_test_heat[3]/Thal_VPLM_L_BF_test_heat[1]
Thal_VPLM_L_BF_no_int_press = Thal_VPLM_L_BF_test_press[3]/Thal_VPLM_L_BF_test_press[1]
Thal_VPLM_R_BF_no_int_heat = Thal_VPLM_R_BF_test_heat[3]/Thal_VPLM_R_BF_test_heat[1]
Thal_VPLM_R_BF_no_int_press = Thal_VPLM_R_BF_test_press[3]/Thal_VPLM_R_BF_test_press[1]
Thal_MD_BF_no_int_heat = Thal_MD_BF_test_heat[3]/Thal_MD_BF_test_heat[1]
Thal_MD_BF_no_int_press = Thal_MD_BF_test_press[3]/Thal_MD_BF_test_press[1]
PAG_BF_no_int_heat = PAG_BF_test_heat[3]/PAG_BF_test_heat[1]
PAG_BF_no_int_press = PAG_BF_test_press[3]/PAG_BF_test_press[1]

list_nociception_bf_no_int_heat = lapply(list(NPS = NPS_BF_no_int_heat,
                                                    dpIns_L = dpIns_L_BF_no_int_heat,      
                                                    dpIns_R = dpIns_R_BF_no_int_heat,
                                                    aMCC_MPFC = aMCC_MPFC_BF_no_int_heat,
                                                    Thal_VPLM_L = Thal_VPLM_L_BF_no_int_heat,
                                                    Thal_VPLM_R = Thal_VPLM_R_BF_no_int_heat,
                                                    Thal_MD = Thal_MD_BF_no_int_heat,
                                                    PAG = PAG_BF_no_int_heat), extractBF)
list_nociception_bf_no_int_press = lapply(list(NPS = NPS_BF_no_int_press,
                                                     dpIns_L = dpIns_L_BF_no_int_press,      
                                                     dpIns_R = dpIns_R_BF_no_int_press,
                                                     aMCC_MPFC = aMCC_MPFC_BF_no_int_press,
                                                     Thal_VPLM_L = Thal_VPLM_L_BF_no_int_press,
                                                     Thal_VPLM_R = Thal_VPLM_R_BF_no_int_press,
                                                     Thal_MD = Thal_MD_BF_no_int_press,
                                                     PAG = PAG_BF_no_int_press), extractBF)

# ultra-wide:
dpIns_L_BF_ultrawide_no_int_heat = dpIns_L_BF_test_ultrawide_heat[3]/dpIns_L_BF_test_ultrawide_heat[1]
dpIns_L_BF_ultrawide_no_int_press = dpIns_L_BF_test_ultrawide_press[3]/dpIns_L_BF_test_ultrawide_press[1]
dpIns_R_BF_ultrawide_no_int_heat = dpIns_R_BF_test_ultrawide_heat[3]/dpIns_R_BF_test_ultrawide_heat[1]
dpIns_R_BF_ultrawide_no_int_press = dpIns_R_BF_test_ultrawide_press[3]/dpIns_R_BF_test_ultrawide_press[1]
aMCC_MPFC_BF_ultrawide_no_int_heat = aMCC_MPFC_BF_test_ultrawide_heat[3]/aMCC_MPFC_BF_test_ultrawide_heat[1]
aMCC_MPFC_BF_ultrawide_no_int_press = aMCC_MPFC_BF_test_ultrawide_press[3]/aMCC_MPFC_BF_test_ultrawide_press[1]
Thal_VPLM_L_BF_ultrawide_no_int_heat = Thal_VPLM_L_BF_test_ultrawide_heat[3]/Thal_VPLM_L_BF_test_ultrawide_heat[1]
Thal_VPLM_L_BF_ultrawide_no_int_press = Thal_VPLM_L_BF_test_ultrawide_press[3]/Thal_VPLM_L_BF_test_ultrawide_press[1]
Thal_VPLM_R_BF_ultrawide_no_int_heat = Thal_VPLM_R_BF_test_ultrawide_heat[3]/Thal_VPLM_R_BF_test_ultrawide_heat[1]
Thal_VPLM_R_BF_ultrawide_no_int_press = Thal_VPLM_R_BF_test_ultrawide_press[3]/Thal_VPLM_R_BF_test_ultrawide_press[1]
Thal_MD_BF_ultrawide_no_int_heat = Thal_MD_BF_test_ultrawide_heat[3]/Thal_MD_BF_test_ultrawide_heat[1]
Thal_MD_BF_ultrawide_no_int_press = Thal_MD_BF_test_ultrawide_press[3]/Thal_MD_BF_test_ultrawide_press[1]
PAG_BF_ultrawide_no_int_heat = PAG_BF_test_ultrawide_heat[3]/PAG_BF_test_ultrawide_heat[1]
PAG_BF_ultrawide_no_int_press = PAG_BF_test_ultrawide_press[3]/PAG_BF_test_ultrawide_press[1]

# put all in a list
list_nociception_bf_ultrawide_no_int_heat = lapply(list(NPS = NPS_BF_ultrawide_no_int_heat,
                                                 dpIns_L = dpIns_L_BF_ultrawide_no_int_heat,
                                                 dpIns_R = dpIns_R_BF_ultrawide_no_int_heat,
                                                 aMCC_MPFC = aMCC_MPFC_BF_ultrawide_no_int_heat,
                                                 Thal_VPLM_L = Thal_VPLM_L_BF_ultrawide_no_int_heat,
                                                 Thal_VPLM_R = Thal_VPLM_R_BF_ultrawide_no_int_heat,
                                                 Thal_MD = Thal_MD_BF_ultrawide_no_int_heat,
                                                 PAG = PAG_BF_ultrawide_no_int_heat), extractBF)
list_nociception_bf_ultrawide_no_int_press = lapply(list(NPS = NPS_BF_ultrawide_no_int_press,
                                                   dpIns_L = dpIns_L_BF_ultrawide_no_int_press,
                                                   dpIns_R = dpIns_R_BF_ultrawide_no_int_press,
                                                   aMCC_MPFC = aMCC_MPFC_BF_ultrawide_no_int_press,
                                                   Thal_VPLM_L = Thal_VPLM_L_BF_ultrawide_no_int_press,
                                                   Thal_VPLM_R = Thal_VPLM_R_BF_ultrawide_no_int_press,
                                                   Thal_MD = Thal_MD_BF_ultrawide_no_int_press,
                                                   PAG = PAG_BF_ultrawide_no_int_press), extractBF)

# organize BF results for the NPS and a-priori nociceptive ROIs in a table
regions = c("NPS", "dpIns_L", "dpIns_R", "aMCC_MPFC", "Thal_VPLM_L", "Thal_VPLM_R", "Thal_MD", "PAG")

heat_wide = array()
press_wide = array()
heat_wide_int = array()
press_wide_int = array()
heat_reg = array()
press_reg = array()
heat_ultrawide = array()
press_ultrawide = array()

for (ind in 1:length(regions)) {
  cur_region = regions[ind]
  heat_wide[ind] = paste(round(list_nociception_bf_wide_no_int_heat[[cur_region]]$bf,3), " (", round(list_nociception_bf_wide_no_int_heat[[cur_region]]$error,3), ")", sep="")
  press_wide[ind] = paste(round(list_nociception_bf_wide_no_int_press[[cur_region]]$bf,3), " (", round(list_nociception_bf_wide_no_int_press[[cur_region]]$error,3), ")", sep="")
  heat_wide_int[ind] = paste(round(list_nociception_bf_wide_int_heat[[cur_region]]$bf,3), " (", round(list_nociception_bf_wide_int_heat[[cur_region]]$error,3), ")", sep="")
  press_wide_int[ind] = paste(round(list_nociception_bf_wide_int_press[[cur_region]]$bf,3), " (", round(list_nociception_bf_wide_int_press[[cur_region]]$error,3), ")", sep="")
  heat_reg[ind] = paste(round(list_nociception_bf_no_int_heat[[cur_region]]$bf,3), " (", round(list_nociception_bf_no_int_heat[[cur_region]]$error,3), ")", sep="")
  press_reg[ind] = paste(round(list_nociception_bf_no_int_press[[cur_region]]$bf,3), " (", round(list_nociception_bf_no_int_press[[cur_region]]$error,3), ")", sep="")
  heat_ultrawide[ind] = paste(round(list_nociception_bf_ultrawide_no_int_heat[[cur_region]]$bf,3), " (", round(list_nociception_bf_ultrawide_no_int_heat[[cur_region]]$error,3), ")", sep="")
  press_ultrawide[ind] = paste(round(list_nociception_bf_ultrawide_no_int_press[[cur_region]]$bf,3), " (", round(list_nociception_bf_ultrawide_no_int_press[[cur_region]]$error,3), ")", sep="")
}

BF_table_nociception = data.frame(region = regions, heat_wide, heat_wide_int, heat_reg, heat_ultrawide, press_wide, press_wide_int, press_reg, press_ultrawide)
write.csv(BF_table_nociception, file = paste("BF_table_nociception", Sys.Date(),".csv", sep = ""))

## higher level regions (expected increase)
#DLPFC (Ctx_p9_46v_R, Ctx_8C_L, Ctx_46_L), mid lateral OFC (Ctx_a10p_R, Ctx_p10p_R), lateral OFC (Ctx_a47r_R), and NAc (V_Striatum_L, V_Striatum_R) 

# dlpfc - left1
dlpfc_L_BF_test_heat = anovaBF(canlab2018_mean_Ctx_8C_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
dlpfc_L_BF_test_press = anovaBF(canlab2018_mean_Ctx_8C_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
dlpfc_L_BF_test_wide_heat = anovaBF(canlab2018_mean_Ctx_8C_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
dlpfc_L_BF_test_wide_press = anovaBF(canlab2018_mean_Ctx_8C_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
dlpfc_L_BF_test_ultrawide_heat = anovaBF(canlab2018_mean_Ctx_8C_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")
dlpfc_L_BF_test_ultrawide_press = anovaBF(canlab2018_mean_Ctx_8C_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")

# dlpfc - right
dlpfc_R1_BF_test_heat = anovaBF(canlab2018_mean_Ctx_p9_46v_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
dlpfc_R1_BF_test_press = anovaBF(canlab2018_mean_Ctx_p9_46v_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
dlpfc_R1_BF_test_wide_heat = anovaBF(canlab2018_mean_Ctx_p9_46v_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
dlpfc_R1_BF_test_wide_press = anovaBF(canlab2018_mean_Ctx_p9_46v_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
dlpfc_R1_BF_test_ultrawide_heat = anovaBF(canlab2018_mean_Ctx_p9_46v_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")
dlpfc_R1_BF_test_ultrawide_press = anovaBF(canlab2018_mean_Ctx_p9_46v_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")

# dlpfc - left2
dlpfc_L2_BF_test_heat = anovaBF(canlab2018_mean_Ctx_46_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
dlpfc_L2_BF_test_press = anovaBF(canlab2018_mean_Ctx_46_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
dlpfc_L2_BF_test_wide_heat = anovaBF(canlab2018_mean_Ctx_46_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
dlpfc_L2_BF_test_wide_press = anovaBF(canlab2018_mean_Ctx_46_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
dlpfc_L2_BF_test_ultrawide_heat = anovaBF(canlab2018_mean_Ctx_46_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")
dlpfc_L2_BF_test_ultrawide_press = anovaBF(canlab2018_mean_Ctx_46_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")

# Anterior OFC 1
anterior_ofc1_BF_test_heat = anovaBF(canlab2018_mean_Ctx_a10p_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
anterior_ofc1_BF_test_press = anovaBF(canlab2018_mean_Ctx_a10p_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
anterior_ofc1_BF_test_wide_heat = anovaBF(canlab2018_mean_Ctx_a10p_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
anterior_ofc1_BF_test_wide_press = anovaBF(canlab2018_mean_Ctx_a10p_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
anterior_ofc1_BF_test_ultrawide_heat = anovaBF(canlab2018_mean_Ctx_a10p_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")
anterior_ofc1_BF_test_ultrawide_press = anovaBF(canlab2018_mean_Ctx_a10p_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")

# Anterior OFC 2
anterior_ofc2_BF_test_heat = anovaBF(canlab2018_mean_Ctx_p10p_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
anterior_ofc2_BF_test_press = anovaBF(canlab2018_mean_Ctx_p10p_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
anterior_ofc2_BF_test_wide_heat = anovaBF(canlab2018_mean_Ctx_p10p_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
anterior_ofc2_BF_test_wide_press = anovaBF(canlab2018_mean_Ctx_p10p_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
anterior_ofc2_BF_test_ultrawide_heat = anovaBF(canlab2018_mean_Ctx_p10p_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")
anterior_ofc2_BF_test_ultrawide_press = anovaBF(canlab2018_mean_Ctx_p10p_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")

# left mid-lateral OFC 1
left_mid_lat_ofc_BF_test_heat = anovaBF(canlab2018_mean_Ctx_11l_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
left_mid_lat_ofc_BF_test_press = anovaBF(canlab2018_mean_Ctx_11l_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
left_mid_lat_ofc_BF_test_wide_heat = anovaBF(canlab2018_mean_Ctx_11l_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
left_mid_lat_ofc_BF_test_wide_press = anovaBF(canlab2018_mean_Ctx_11l_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
left_mid_lat_ofc_BF_test_ultrawide_heat = anovaBF(canlab2018_mean_Ctx_11l_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")
left_mid_lat_ofc_BF_test_ultrawide_press = anovaBF(canlab2018_mean_Ctx_11l_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")

# mid-lateral OFC 2
right_mid_lat_ofc_BF_test_heat = anovaBF(canlab2018_mean_Ctx_11l_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
right_mid_lat_ofc_BF_test_press = anovaBF(canlab2018_mean_Ctx_11l_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
right_mid_lat_ofc_BF_test_wide_heat = anovaBF(canlab2018_mean_Ctx_11l_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
right_mid_lat_ofc_BF_test_wide_press = anovaBF(canlab2018_mean_Ctx_11l_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
right_mid_lat_ofc_BF_test_ultrawide_heat = anovaBF(canlab2018_mean_Ctx_11l_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")
right_mid_lat_ofc_BF_test_ultrawide_press = anovaBF(canlab2018_mean_Ctx_11l_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")

# lateral OFC
lateral_ofc_BF_test_heat = anovaBF(canlab2018_mean_Ctx_a47r_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
lateral_ofc_BF_test_press = anovaBF(canlab2018_mean_Ctx_a47r_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
lateral_ofc_BF_test_wide_heat = anovaBF(canlab2018_mean_Ctx_a47r_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
lateral_ofc_BF_test_wide_press = anovaBF(canlab2018_mean_Ctx_a47r_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
lateral_ofc_BF_test_ultrawide_heat = anovaBF(canlab2018_mean_Ctx_a47r_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")
lateral_ofc_BF_test_ultrawide_press = anovaBF(canlab2018_mean_Ctx_a47r_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")

# NAc - right
nac_R_BF_test_heat = anovaBF(canlab2018_mean_V_Striatum_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
nac_R_BF_test_press = anovaBF(canlab2018_mean_V_Striatum_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
nac_R_BF_test_wide_heat = anovaBF(canlab2018_mean_V_Striatum_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
nac_R_BF_test_wide_press = anovaBF(canlab2018_mean_V_Striatum_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
nac_R_BF_test_ultrawide_heat = anovaBF(canlab2018_mean_V_Striatum_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")
nac_R_BF_test_ultrawide_press = anovaBF(canlab2018_mean_V_Striatum_R ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")

# NAc - left
nac_L_BF_test_heat = anovaBF(canlab2018_mean_V_Striatum_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
nac_L_BF_test_press = anovaBF(canlab2018_mean_V_Striatum_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain")
nac_L_BF_test_wide_heat = anovaBF(canlab2018_mean_V_Striatum_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
nac_L_BF_test_wide_press = anovaBF(canlab2018_mean_V_Striatum_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "wide")
nac_L_BF_test_ultrawide_heat = anovaBF(canlab2018_mean_V_Striatum_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")
nac_L_BF_test_ultrawide_press = anovaBF(canlab2018_mean_V_Striatum_L ~ stimLvl + prodicaine + participant_ID + family_ID, data = data[data$heat==-0.5,], whichRandom = c("participant_ID", "family_ID"), whichModels="withmain", rscaleFixed = "ultrawide")

# BFs for a-priori higher-level ROIs, wide, no interaction
dlpfc_L_BF_wide_no_int_heat = dlpfc_L_BF_test_wide_heat[3]/dlpfc_L_BF_test_wide_heat[1]
dlpfc_L_BF_wide_no_int_press = dlpfc_L_BF_test_wide_press[3]/dlpfc_L_BF_test_wide_press[1]
dlpfc_L2_BF_wide_no_int_heat = dlpfc_L2_BF_test_wide_heat[3]/dlpfc_L2_BF_test_wide_heat[1]
dlpfc_L2_BF_wide_no_int_press = dlpfc_L2_BF_test_wide_press[3]/dlpfc_L2_BF_test_wide_press[1]
dlpfc_R1_BF_wide_no_int_heat = dlpfc_R1_BF_test_wide_heat[3]/dlpfc_R1_BF_test_wide_heat[1]
dlpfc_R1_BF_wide_no_int_press = dlpfc_R1_BF_test_wide_press[3]/dlpfc_R1_BF_test_wide_press[1]
lateral_ofc_BF_wide_no_int_heat = lateral_ofc_BF_test_wide_heat[3]/lateral_ofc_BF_test_wide_heat[1]
lateral_ofc_BF_wide_no_int_press = lateral_ofc_BF_test_wide_press[3]/lateral_ofc_BF_test_wide_press[1]
anterior_ofc1_BF_wide_no_int_heat = anterior_ofc1_BF_test_wide_heat[3]/anterior_ofc1_BF_test_wide_heat[1]
anterior_ofc1_BF_wide_no_int_press = anterior_ofc1_BF_test_wide_press[3]/anterior_ofc1_BF_test_wide_press[1]
anterior_ofc2_BF_wide_no_int_heat = anterior_ofc2_BF_test_wide_heat[3]/anterior_ofc2_BF_test_wide_heat[1]
anterior_ofc2_BF_wide_no_int_press = anterior_ofc2_BF_test_wide_press[3]/anterior_ofc2_BF_test_wide_press[1]
nac_R_BF_wide_no_int_heat = nac_R_BF_test_wide_heat[3]/nac_R_BF_test_wide_heat[1]
nac_R_BF_wide_no_int_press = nac_R_BF_test_wide_press[3]/nac_R_BF_test_wide_press[1]
nac_L_BF_wide_no_int_heat = nac_L_BF_test_wide_heat[3]/nac_L_BF_test_wide_heat[1]
nac_L_BF_wide_no_int_press = nac_L_BF_test_wide_press[3]/nac_L_BF_test_wide_press[1]
right_mid_lat_ofc_BF_wide_no_int_heat = right_mid_lat_ofc_BF_test_wide_heat[3]/right_mid_lat_ofc_BF_test_wide_heat[1]
right_mid_lat_ofc_BF_wide_no_int_press = right_mid_lat_ofc_BF_test_wide_press[3]/right_mid_lat_ofc_BF_test_wide_press[1]
left_mid_lat_ofc_BF_wide_no_int_heat = left_mid_lat_ofc_BF_test_wide_heat[3]/left_mid_lat_ofc_BF_test_wide_heat[1]
left_mid_lat_ofc_BF_wide_no_int_press = left_mid_lat_ofc_BF_test_wide_press[3]/left_mid_lat_ofc_BF_test_wide_press[1]

# put all in a list
list_higher_processing_bf_wide_no_int_heat = lapply(list(SIIPS = siips_BF_wide_no_int_heat,
                                             dlpfc_L = dlpfc_L_BF_wide_no_int_heat,
                                             dlpfc_L2 = dlpfc_L2_BF_wide_no_int_heat,
                                             dlpfc_R1 = dlpfc_R1_BF_wide_no_int_heat,
                                             lateral_ofc = lateral_ofc_BF_wide_no_int_heat,
                                             anterior_ofc1 = anterior_ofc1_BF_wide_no_int_heat,
                                             anterior_ofc2 = anterior_ofc2_BF_wide_no_int_heat,
                                             nac_R = nac_R_BF_wide_no_int_heat,
                                             nac_L = nac_L_BF_wide_no_int_heat,
                                             right_mid_lat_ofc = right_mid_lat_ofc_BF_wide_no_int_heat,
                                             left_mid_lat_ofc = left_mid_lat_ofc_BF_wide_no_int_heat), extractBF)
list_higher_processing_bf_wide_no_int_press = lapply(list(SIIPS = siips_BF_wide_no_int_press,
                                             dlpfc_L = dlpfc_L_BF_wide_no_int_press,
                                             dlpfc_L2 = dlpfc_L2_BF_wide_no_int_press,
                                             dlpfc_R1 = dlpfc_R1_BF_wide_no_int_press,
                                             lateral_ofc = lateral_ofc_BF_wide_no_int_press,
                                             anterior_ofc1 = anterior_ofc1_BF_wide_no_int_press,
                                             anterior_ofc2 = anterior_ofc2_BF_wide_no_int_press,
                                             nac_R = nac_R_BF_wide_no_int_press,
                                             nac_L = nac_L_BF_wide_no_int_press,
                                             right_mid_lat_ofc = right_mid_lat_ofc_BF_wide_no_int_press,
                                             left_mid_lat_ofc = left_mid_lat_ofc_BF_wide_no_int_press), extractBF)

# robustness tests:
# with interaction
dlpfc_L_BF_wide_int_heat = dlpfc_L_BF_test_wide_heat[4]/dlpfc_L_BF_test_wide_heat[1]
dlpfc_L_BF_wide_int_press = dlpfc_L_BF_test_wide_press[4]/dlpfc_L_BF_test_wide_press[1]
dlpfc_L2_BF_wide_int_heat = dlpfc_L2_BF_test_wide_heat[4]/dlpfc_L2_BF_test_wide_heat[1]
dlpfc_L2_BF_wide_int_press = dlpfc_L2_BF_test_wide_press[4]/dlpfc_L2_BF_test_wide_press[1]
dlpfc_R1_BF_wide_int_heat = dlpfc_R1_BF_test_wide_heat[4]/dlpfc_R1_BF_test_wide_heat[1]
dlpfc_R1_BF_wide_int_press = dlpfc_R1_BF_test_wide_press[4]/dlpfc_R1_BF_test_wide_press[1]
lateral_ofc_BF_wide_int_heat = lateral_ofc_BF_test_wide_heat[4]/lateral_ofc_BF_test_wide_heat[1]
lateral_ofc_BF_wide_int_press = lateral_ofc_BF_test_wide_press[4]/lateral_ofc_BF_test_wide_press[1]
anterior_ofc1_BF_wide_int_heat = anterior_ofc1_BF_test_wide_heat[4]/anterior_ofc1_BF_test_wide_heat[1]
anterior_ofc1_BF_wide_int_press = anterior_ofc1_BF_test_wide_press[4]/anterior_ofc1_BF_test_wide_press[1]
anterior_ofc2_BF_wide_int_heat = anterior_ofc2_BF_test_wide_heat[4]/anterior_ofc2_BF_test_wide_heat[1]
anterior_ofc2_BF_wide_int_press = anterior_ofc2_BF_test_wide_press[4]/anterior_ofc2_BF_test_wide_press[1]
nac_R_BF_wide_int_heat = nac_R_BF_test_wide_heat[4]/nac_R_BF_test_wide_heat[1]
nac_R_BF_wide_int_press = nac_R_BF_test_wide_press[4]/nac_R_BF_test_wide_press[1]
nac_L_BF_wide_int_heat = nac_L_BF_test_wide_heat[4]/nac_L_BF_test_wide_heat[1]
nac_L_BF_wide_int_press = nac_L_BF_test_wide_press[4]/nac_L_BF_test_wide_press[1]
right_mid_lat_ofc_BF_wide_int_heat = right_mid_lat_ofc_BF_test_wide_heat[4]/right_mid_lat_ofc_BF_test_wide_heat[1]
right_mid_lat_ofc_BF_wide_int_press = right_mid_lat_ofc_BF_test_wide_press[4]/right_mid_lat_ofc_BF_test_wide_press[1]
left_mid_lat_ofc_BF_wide_int_heat = left_mid_lat_ofc_BF_test_wide_heat[4]/left_mid_lat_ofc_BF_test_wide_heat[1]
left_mid_lat_ofc_BF_wide_int_press = left_mid_lat_ofc_BF_test_wide_press[4]/left_mid_lat_ofc_BF_test_wide_press[1]

# put all in a list
list_higher_processing_bf_wide_int_heat = lapply(list(SIIPS = siips_BF_wide_int_heat,
                                             dlpfc_L = dlpfc_L_BF_wide_int_heat,
                                             dlpfc_L2 = dlpfc_L2_BF_wide_int_heat,
                                             dlpfc_R1 = dlpfc_R1_BF_wide_int_heat,
                                             lateral_ofc = lateral_ofc_BF_wide_int_heat,
                                             anterior_ofc1 = anterior_ofc1_BF_wide_int_heat,
                                             anterior_ofc2 = anterior_ofc2_BF_wide_int_heat,
                                             nac_R = nac_R_BF_wide_int_heat,
                                             nac_L = nac_L_BF_wide_int_heat,
                                             left_mid_lat_ofc = left_mid_lat_ofc_BF_wide_int_heat,
                                             right_mid_lat_ofc = right_mid_lat_ofc_BF_wide_int_heat), extractBF)
list_higher_processing_bf_wide_int_press = lapply(list(SIIPS = siips_BF_wide_int_press,
                                              dlpfc_L = dlpfc_L_BF_wide_int_press,
                                              dlpfc_L2 = dlpfc_L2_BF_wide_int_press,
                                              dlpfc_R1 = dlpfc_R1_BF_wide_int_press,
                                              lateral_ofc = lateral_ofc_BF_wide_int_press,
                                              anterior_ofc1 = anterior_ofc1_BF_wide_int_press,
                                              anterior_ofc2 = anterior_ofc2_BF_wide_int_press,
                                              nac_R = nac_R_BF_wide_int_press,
                                              nac_L = nac_L_BF_wide_int_press,
                                              left_mid_lat_ofc = left_mid_lat_ofc_BF_wide_int_press,
                                              right_mid_lat_ofc = right_mid_lat_ofc_BF_wide_int_press), extractBF)

# default width:
dlpfc_L_BF_no_int_heat = dlpfc_L_BF_test_heat[3]/dlpfc_L_BF_test_heat[1]
dlpfc_L_BF_no_int_press = dlpfc_L_BF_test_press[3]/dlpfc_L_BF_test_press[1]
dlpfc_L2_BF_no_int_heat = dlpfc_L2_BF_test_heat[3]/dlpfc_L2_BF_test_heat[1]
dlpfc_L2_BF_no_int_press = dlpfc_L2_BF_test_press[3]/dlpfc_L2_BF_test_press[1]
dlpfc_R1_BF_no_int_heat = dlpfc_R1_BF_test_heat[3]/dlpfc_R1_BF_test_heat[1]
dlpfc_R1_BF_no_int_press = dlpfc_R1_BF_test_press[3]/dlpfc_R1_BF_test_press[1]
lateral_ofc_BF_no_int_heat = lateral_ofc_BF_test_heat[3]/lateral_ofc_BF_test_heat[1]
lateral_ofc_BF_no_int_press = lateral_ofc_BF_test_press[3]/lateral_ofc_BF_test_press[1]
anterior_ofc1_BF_no_int_heat = anterior_ofc1_BF_test_heat[3]/anterior_ofc1_BF_test_heat[1]
anterior_ofc1_BF_no_int_press = anterior_ofc1_BF_test_press[3]/anterior_ofc1_BF_test_press[1]
anterior_ofc2_BF_no_int_heat = anterior_ofc2_BF_test_heat[3]/anterior_ofc2_BF_test_heat[1]
anterior_ofc2_BF_no_int_press = anterior_ofc2_BF_test_press[3]/anterior_ofc2_BF_test_press[1]
nac_R_BF_no_int_heat = nac_R_BF_test_heat[3]/nac_R_BF_test_heat[1]
nac_R_BF_no_int_press = nac_R_BF_test_press[3]/nac_R_BF_test_press[1]
nac_L_BF_no_int_heat = nac_L_BF_test_heat[3]/nac_L_BF_test_heat[1]
nac_L_BF_no_int_press = nac_L_BF_test_press[3]/nac_L_BF_test_press[1]
right_mid_lat_ofc_BF_no_int_heat = right_mid_lat_ofc_BF_test_heat[3]/right_mid_lat_ofc_BF_test_heat[1]
right_mid_lat_ofc_BF_no_int_press = right_mid_lat_ofc_BF_test_press[3]/right_mid_lat_ofc_BF_test_press[1]
left_mid_lat_ofc_BF_no_int_heat = left_mid_lat_ofc_BF_test_heat[3]/left_mid_lat_ofc_BF_test_heat[1]
left_mid_lat_ofc_BF_no_int_press = left_mid_lat_ofc_BF_test_press[3]/left_mid_lat_ofc_BF_test_press[1]

list_higher_processing_bf_no_int_heat = lapply(list(SIIPS = siips_BF_no_int_heat,
                                             dlpfc_L = dlpfc_L_BF_no_int_heat,
                                             dlpfc_L2 = dlpfc_L2_BF_no_int_heat,
                                             dlpfc_R1 = dlpfc_R1_BF_no_int_heat,
                                             lateral_ofc = lateral_ofc_BF_no_int_heat,
                                             anterior_ofc1 = anterior_ofc1_BF_no_int_heat,
                                             anterior_ofc2 = anterior_ofc2_BF_no_int_heat,
                                             nac_R = nac_R_BF_no_int_heat,
                                             nac_L = nac_L_BF_no_int_heat,
                                             left_mid_lat_ofc = left_mid_lat_ofc_BF_no_int_heat,
                                             right_mid_lat_ofc = right_mid_lat_ofc_BF_no_int_heat), extractBF)
list_higher_processing_bf_no_int_press = lapply(list(SIIPS = siips_BF_no_int_press,
                                              dlpfc_L = dlpfc_L_BF_no_int_press,
                                              dlpfc_L2 = dlpfc_L2_BF_no_int_press,
                                              dlpfc_R1 = dlpfc_R1_BF_no_int_press,
                                              lateral_ofc = lateral_ofc_BF_no_int_press,
                                              anterior_ofc1 = anterior_ofc1_BF_no_int_press,
                                              anterior_ofc2 = anterior_ofc2_BF_no_int_press,
                                              nac_R = nac_R_BF_no_int_press,
                                              nac_L = nac_L_BF_no_int_press,
                                              left_mid_lat_ofc = left_mid_lat_ofc_BF_no_int_press,
                                              right_mid_lat_ofc = right_mid_lat_ofc_BF_no_int_press), extractBF)

# ultra-wide:
dlpfc_L_BF_ultrawide_no_int_heat = dlpfc_L_BF_test_ultrawide_heat[3]/dlpfc_L_BF_test_ultrawide_heat[1]
dlpfc_L_BF_ultrawide_no_int_press = dlpfc_L_BF_test_ultrawide_press[3]/dlpfc_L_BF_test_ultrawide_press[1]
dlpfc_L2_BF_ultrawide_no_int_heat = dlpfc_L2_BF_test_ultrawide_heat[3]/dlpfc_L2_BF_test_ultrawide_heat[1]
dlpfc_L2_BF_ultrawide_no_int_press = dlpfc_L2_BF_test_ultrawide_press[3]/dlpfc_L2_BF_test_ultrawide_press[1]
dlpfc_R1_BF_ultrawide_no_int_heat = dlpfc_R1_BF_test_ultrawide_heat[3]/dlpfc_R1_BF_test_ultrawide_heat[1]
dlpfc_R1_BF_ultrawide_no_int_press = dlpfc_R1_BF_test_ultrawide_press[3]/dlpfc_R1_BF_test_ultrawide_press[1]
lateral_ofc_BF_ultrawide_no_int_heat = lateral_ofc_BF_test_ultrawide_heat[3]/lateral_ofc_BF_test_ultrawide_heat[1]
lateral_ofc_BF_ultrawide_no_int_press = lateral_ofc_BF_test_ultrawide_press[3]/lateral_ofc_BF_test_ultrawide_press[1]
anterior_ofc1_BF_ultrawide_no_int_heat = anterior_ofc1_BF_test_ultrawide_heat[3]/anterior_ofc1_BF_test_ultrawide_heat[1]
anterior_ofc1_BF_ultrawide_no_int_press = anterior_ofc1_BF_test_ultrawide_press[3]/anterior_ofc1_BF_test_ultrawide_press[1]
anterior_ofc2_BF_ultrawide_no_int_heat = anterior_ofc2_BF_test_ultrawide_heat[3]/anterior_ofc2_BF_test_ultrawide_heat[1]
anterior_ofc2_BF_ultrawide_no_int_press = anterior_ofc2_BF_test_ultrawide_press[3]/anterior_ofc2_BF_test_ultrawide_press[1]
nac_R_BF_ultrawide_no_int_heat = nac_R_BF_test_ultrawide_heat[3]/nac_R_BF_test_ultrawide_heat[1]
nac_R_BF_ultrawide_no_int_press = nac_R_BF_test_ultrawide_press[3]/nac_R_BF_test_ultrawide_press[1]
nac_L_BF_ultrawide_no_int_heat = nac_L_BF_test_ultrawide_heat[3]/nac_L_BF_test_ultrawide_heat[1]
nac_L_BF_ultrawide_no_int_press = nac_L_BF_test_ultrawide_press[3]/nac_L_BF_test_ultrawide_press[1]
right_mid_lat_ofc_BF_ultrawide_no_int_heat = right_mid_lat_ofc_BF_test_ultrawide_heat[3]/right_mid_lat_ofc_BF_test_ultrawide_heat[1]
right_mid_lat_ofc_BF_ultrawide_no_int_press = right_mid_lat_ofc_BF_test_ultrawide_press[3]/right_mid_lat_ofc_BF_test_ultrawide_press[1]
left_mid_lat_ofc_BF_ultrawide_no_int_heat = left_mid_lat_ofc_BF_test_ultrawide_heat[3]/left_mid_lat_ofc_BF_test_ultrawide_heat[1]
left_mid_lat_ofc_BF_ultrawide_no_int_press = left_mid_lat_ofc_BF_test_ultrawide_press[3]/left_mid_lat_ofc_BF_test_ultrawide_press[1]

list_higher_processing_bf_ultrawide_no_int_heat = lapply(list(SIIPS = siips_BF_ultrawide_no_int_heat,
                                             dlpfc_L = dlpfc_L_BF_ultrawide_no_int_heat,
                                             dlpfc_L2 = dlpfc_L2_BF_ultrawide_no_int_heat,
                                             dlpfc_R1 = dlpfc_R1_BF_ultrawide_no_int_heat,
                                             lateral_ofc = lateral_ofc_BF_ultrawide_no_int_heat,
                                             anterior_ofc1 = anterior_ofc1_BF_ultrawide_no_int_heat,
                                             anterior_ofc2 = anterior_ofc2_BF_ultrawide_no_int_heat,
                                             nac_R = nac_R_BF_ultrawide_no_int_heat,
                                             nac_L = nac_L_BF_ultrawide_no_int_heat,
                                             left_mid_lat_ofc = left_mid_lat_ofc_BF_ultrawide_no_int_heat,
                                             right_mid_lat_ofc = right_mid_lat_ofc_BF_ultrawide_no_int_heat), extractBF)
list_higher_processing_bf_ultrawide_no_int_press = lapply(list(SIIPS = siips_BF_ultrawide_no_int_press,
                                              dlpfc_L = dlpfc_L_BF_ultrawide_no_int_press,
                                              dlpfc_L2 = dlpfc_L2_BF_ultrawide_no_int_press,
                                              dlpfc_R1 = dlpfc_R1_BF_ultrawide_no_int_press,
                                              lateral_ofc = lateral_ofc_BF_ultrawide_no_int_press,
                                              anterior_ofc1 = anterior_ofc1_BF_ultrawide_no_int_press,
                                              anterior_ofc2 = anterior_ofc2_BF_ultrawide_no_int_press,
                                              nac_R = nac_R_BF_ultrawide_no_int_press,
                                              nac_L = nac_L_BF_ultrawide_no_int_press,
                                              left_mid_lat_ofc = left_mid_lat_ofc_BF_ultrawide_no_int_press,
                                              right_mid_lat_ofc = right_mid_lat_ofc_BF_ultrawide_no_int_press), extractBF)

# organize BF results for a-priori ROIs (expected increase) in a table
regions = c("SIIPS", "dlpfc_L", "dlpfc_L2", "dlpfc_R1", "lateral_ofc", "anterior_ofc1", "anterior_ofc2", "nac_L", "nac_R", "left_mid_lat_ofc", "right_mid_lat_ofc")

heat_wide = array()
press_wide = array()
heat_wide_int = array()
press_wide_int = array()
heat_reg = array()
press_reg = array()
heat_ultrawide = array()
press_ultrawide = array()

for (ind in 1:length(regions)) {
  cur_region = regions[ind]
  heat_wide[ind] = paste(round(list_higher_processing_bf_wide_no_int_heat[[cur_region]]$bf,3), " (", round(list_higher_processing_bf_wide_no_int_heat[[cur_region]]$error,3), ")", sep="")
  press_wide[ind] = paste(round(list_higher_processing_bf_wide_no_int_press[[cur_region]]$bf,3), " (", round(list_higher_processing_bf_wide_no_int_press[[cur_region]]$error,3), ")", sep="")
  heat_wide_int[ind] = paste(round(list_higher_processing_bf_wide_int_heat[[cur_region]]$bf,3), " (", round(list_higher_processing_bf_wide_int_heat[[cur_region]]$error,3), ")", sep="")
  press_wide_int[ind] = paste(round(list_higher_processing_bf_wide_int_press[[cur_region]]$bf,3), " (", round(list_higher_processing_bf_wide_int_press[[cur_region]]$error,3), ")", sep="")
  heat_reg[ind] = paste(round(list_higher_processing_bf_no_int_heat[[cur_region]]$bf,3), " (", round(list_higher_processing_bf_no_int_heat[[cur_region]]$error,3), ")", sep="")
  press_reg[ind] = paste(round(list_higher_processing_bf_no_int_press[[cur_region]]$bf,3), " (", round(list_higher_processing_bf_no_int_press[[cur_region]]$error,3), ")", sep="")
  heat_ultrawide[ind] = paste(round(list_higher_processing_bf_ultrawide_no_int_heat[[cur_region]]$bf,3), " (", round(list_higher_processing_bf_ultrawide_no_int_heat[[cur_region]]$error,3), ")", sep="")
  press_ultrawide[ind] = paste(round(list_higher_processing_bf_ultrawide_no_int_press[[cur_region]]$bf,3), " (", round(list_higher_processing_bf_ultrawide_no_int_press[[cur_region]]$error,3), ")", sep="")
}

BF_table_higher_processing = data.frame(region = regions, heat_wide, heat_wide_int, heat_reg, heat_ultrawide, press_wide, press_wide_int, press_reg, press_ultrawide)
write.csv(BF_table_higher_processing, file = paste("BF_table_higher_processing", Sys.Date(),".csv", sep = ""))

# save workspace
save.image(file = paste("BF_tests_", Sys.Date(), ".RData", sep = ""))