library(lme4)
library(dplyr)
library(bayestestR)

fullTrain2 <- read.csv("Study 2/Cleaning/output/fullTrain.csv")
indDiffs2 <- read.csv("Study 2/Cleaning/output/fullTest.csv") %>% filter(!duplicated(subID))

# Merge indDiffs2 to fullTrain2
# The variables are likely inOutTherm (ingroup bias) and MGIS (social identification)
if (!"inOutTherm" %in% names(indDiffs2)) {
  # Maybe it's calculated as InUCLATherm/InCSULATherm or Therm_1 - Therm_2?
  # Let's check available columns
  print("inOutTherm not found in indDiffs2")
} else {
  fullTrain2 <- fullTrain2 %>% left_join(indDiffs2 %>% select(subID, groupHomoph, inOutTherm, MGIS), by="subID")
}

fullTrain2$inGsimS.Z <- scale(fullTrain2$inGsimS)
fullTrain2$outGsimS.Z <- scale(fullTrain2$outGsimS)

# 1. Similarity-to-outgroup and ingroup bias (inOutTherm)
if ("inOutTherm" %in% names(fullTrain2)) {
  fullTrain2$inOutTherm.Z <- scale(fullTrain2$inOutTherm)
  m1 <- lmer(scale(selfResp) ~ inGsimS.Z * inOutTherm.Z + outGsimS.Z * inOutTherm.Z + (1 | subID) + (1 | trait), data=fullTrain2)
  print("Model 1: Ingroup Bias Interaction")
  print(summary(m1)$coefficients)
}

# 2. Ingroup similarity and MGIS by outgroup
if ("MGIS" %in% names(fullTrain2)) {
  fullTrain2$MGIS.Z <- scale(fullTrain2$MGIS)
  fullTrain2$outgroup <- factor(fullTrain2$outgroup)
  m2 <- lmer(scale(selfResp) ~ inGsimS.Z * outgroup * MGIS.Z + (1 | subID) + (1 | trait), data=fullTrain2)
  print("Model 2: MGIS and outgroup Interaction")
  print(summary(m2)$coefficients)
}

# 3. Segregation (groupHomoph) by outgroup
indDiffs2$groupHomoph.Z <- scale(indDiffs2$groupHomoph)
indDiffs2$outgroup <- factor(indDiffs2$outgroup)
m3 <- lm(groupHomoph.Z ~ outgroup, data=indDiffs2)
print("Model 3: groupHomoph by outgroup")
print(summary(m3)$coefficients)
