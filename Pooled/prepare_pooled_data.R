# Prepare Pooled Data for Stan
library(tidyverse)
library(here)
library(igraph)

# 1. Load Study 1
message("Loading Study 1...")
s1_test <- read.csv("Study 1/Cleaning/output/fullTest.csv") %>% filter(!is.na(ingChoiceN))
s1_train <- read.csv("Study 1/Cleaning/output/fullTrain.csv") %>% filter(!is.na(selfResp))
common_s1 <- intersect(unique(s1_test$subID), unique(s1_train$subID))
s1_test <- s1_test %>% filter(subID %in% common_s1)
s1_train <- s1_train %>% filter(subID %in% common_s1)
message("S1 Test rows: ", nrow(s1_test), " Train rows: ", nrow(s1_train))

# 2. Load Study 2
message("Loading Study 2...")
s2_test <- read.csv("Study 2/Cleaning/output/fullTest.csv") %>% filter(!is.na(ingChoiceN))
s2_train <- read.csv("Study 2/Cleaning/output/fullTrain_fixed.csv") %>% filter(!is.na(selfResp))
common_s2 <- intersect(unique(s2_test$subID), unique(s2_train$subID))
s2_test <- s2_test %>% filter(subID %in% common_s2)
s2_train <- s2_train %>% filter(subID %in% common_s2)
message("S2 Test rows: ", nrow(s2_test), " Train rows: ", nrow(s2_train))

# 3. Load Study 3
message("Loading Study 3...")
s3_test <- read.csv("Study 3/Cleaning/output/fullTest_fixed.csv") %>% filter(!is.na(ingChoiceN))
s3_train <- read.csv("Study 3/Cleaning/output/fullTrain_fixed.csv") %>% filter(!is.na(selfResp))
common_s3 <- intersect(unique(s3_test$subID), unique(s3_train$subID))
s3_test <- s3_test %>% filter(subID %in% common_s3)
s3_train <- s3_train %>% filter(subID %in% common_s3)
message("S3 Test rows: ", nrow(s3_test), " Train rows: ", nrow(s3_train))

# 4. Combine
s1_test$study <- 1
s2_test$study <- 2
s3_test$study <- 3
s1_train$study <- 1
s2_train$study <- 2
s3_train$study <- 3

fulldf <- bind_rows(
  s1_test %>% select(subID, study, trait, ingChoiceN, Idx, trialTotalT2),
  s2_test %>% select(subID, study, trait, ingChoiceN, Idx, trialTotalT2),
  s3_test %>% select(subID, study, trait, ingChoiceN, Idx, trialTotalT2)
)

traindf <- bind_rows(
  s1_train %>% select(subID, study, trait, selfResp, Idx),
  s2_train %>% select(subID, study, trait, selfResp, Idx),
  s3_train %>% select(subID, study, trait, selfResp, Idx)
)

# Create unique subID per study
fulldf <- fulldf %>% mutate(uniqueSubID = paste0(study, "_", subID))
traindf <- traindf %>% mutate(uniqueSubID = paste0(study, "_", subID))

uIds <- sort(unique(fulldf$uniqueSubID))
maxSubjs <- length(uIds)
maxTrials <- max(fulldf$trialTotalT2, na.rm=TRUE)
maxTrain <- 91

# Adjacency Matrix
posDf <- read.csv("Pooled/input/adjacencyMatrix_p.csv")
simMat <- similarity(graph_from_adjacency_matrix(as.matrix(posDf), mode = "undirected"), method = "dice")

# Prepare Stan Data
stan_data <- list(
  nSubjects = maxSubjs,
  nStudies = 3,
  subjStudy = integer(maxSubjs),
  maxTrials = maxTrials,
  maxTrain = maxTrain,
  nTrain = integer(maxSubjs),
  nTrials = integer(maxSubjs),
  groupChoice = array(0, c(maxSubjs, maxTrials)),
  prevSim = array(0, c(maxSubjs, maxTrials, maxTrain)),
  prevSelf = array(0, c(maxSubjs, maxTrain))
)

message("Populating Stan data...")
for(i in 1:maxSubjs) {
  s_df <- filter(fulldf, uniqueSubID == uIds[i])
  s_train <- filter(traindf, uniqueSubID == uIds[i])
  
  stan_data$subjStudy[i] <- s_df$study[1]
  stan_data$nTrials[i] <- nrow(s_df)
  stan_data$nTrain[i] <- nrow(s_train)
  
  if (nrow(s_df) > 0) {
    stan_data$groupChoice[i, 1:nrow(s_df)] <- s_df$ingChoiceN + 1
    if (nrow(s_train) > 0) {
      stan_data$prevSim[i, 1:nrow(s_df), 1:nrow(s_train)] <- simMat[s_df$Idx, s_train$Idx]
    }
  }
  if (nrow(s_train) > 0) {
    stan_data$prevSelf[i, 1:nrow(s_train)] <- s_train$selfResp
  }
}

saveRDS(stan_data, "pooled_stan_data.rds")
message("Pooled Stan data ready. Total subjects: ", maxSubjs)
